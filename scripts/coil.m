#!/usr/bin/env octave

clear all
close all 
clc

global inductance = 0.0;
global resistance = 0.0;

%% Arguments
LayerWinding = str2num(argv(){1});      % Turns per layer
Layer = str2num(argv(){2});   		% Number of layers
dcu  = str2num(argv(){3});         	% Pure Cu-diameter
current = str2num(argv(){4}); 

strands = 1;
Cduct = 51; %58!!!!!!!
dens = 8.96; %g/cm3
turns = LayerWinding * Layer; % Number of single conductor wires
dx = dcu*1.15;       % distance center conductor to center conductor

%   |A|<-B->|A|                  kana�
%   +-+     +-+     -           |<-B->|
%   | |     | |     |        -  +-----+
%   | |     | |     |        |  |     |
%   | |     | |     C        D  |     |
%   | |     | |     |        |  |     |
%   | |     | |     |        -  +-----+
%   +-+     +-+     -

A = Layer*dx;
B = 48; %40
C = LayerWinding*strands*dx;
D = 145; %129
Depth = D+2*A;

MLT = 2*(A+B)+2*Depth; %mm
Scu = strands*3.14*(dcu/2)^2; %mm2
Wcu = dens * MLT * turns * Scu * 10^-6; %kg



function initFEMM(Depth)
	path(path,'~/.wine/drive_c/Software/femm42/mfiles')   % sets the path
	openfemm                        % launches FEMM
	newdocument(0)                  % create new FEMM M-static file

	mi_probdef(0,  'millimeters','planar', 1e-8, Depth, 28);
	%-- 1e-8=Max.Fehler f�r Solver (Max. Error for Solver)
	%-- 30� =Min.Winkel f�r Netz (max Angle of mesh)
endfunction

function addMaterials(Cduct, strands, dcu)
	wireName = sprintf('%0.1fmm',dcu);
	% mi_addmaterial( �name�,mux,muy, Hc ,   J,Cduct,LamD,PhiHmax,LamFill,LamType,PhiHx,PhiHy,   nStr,dWire);
	mi_addmaterial(   'Air',  1,  1,   0,   0,    0,   0,      0,      0,      0,    0,    0,       0,   0);
	mi_addmaterial(wireName,  1,  1,   0,   0,Cduct,   0,      0,      0,      5,    0,    0, strands, dcu);
endfunction

function setupWorkspace(A, B, C) 
	%--    ---upper---
	%--   |           |
	%--   L           R
	%--   |           |
	%--    ---lower---

	%-- Workspace coordinates
	scale_factor = 8;
	upper = C+B/2*scale_factor;
	lower = -B/2*scale_factor;
	L = -A-B/2-B/2*scale_factor;
	R = -L;
	mi_addnode(L,lower);
	mi_addnode(L,upper);
	mi_addnode(R,upper);
	mi_addnode(R,lower);
	mi_addsegment(L,lower,L,upper);
	mi_addsegment(L,upper,R,upper);
	mi_addsegment(R,upper,R,lower);
	mi_addsegment(R,lower,L,lower);

	%--mi_addboundprop("propname", A0, A1, A2, Phi, Mu, Sig, c0, c1, BdryFormat)
	u0 = 12.57E-7;
	c0 = 1.0 / ( C/1000*u0);
	c1=0;  a0=0;  a1=0; a2=0; Phi=0; Mu=0; Sig=0;
	Bndry = 2;  %Mixed
	mi_addboundprop('NULL', a0, a1, a2, Phi, Mu, Sig, c0, c1, Bndry);

	mi_selectsegment(0,lower);
	mi_selectsegment(0,upper);
	mi_selectsegment(L,0);
	mi_selectsegment(R,0);
	mi_setsegmentprop('NULL', 1, 1, 0, 99);
	mi_clearselected();

	%-- Air Properties
	mi_addblocklabel(lower+10,L+50);
	mi_selectlabel(lower+10,L+50);
	%-- mi_setblockprop("name",automesh,meshsize,"incircuit",magdir,group,turns);
	mi_setblockprop('Air',    1,       5 ,     '1',          0,     1,    0);
	mi_clearselected();

	mi_zoomnatural();
endfunction

function drawCoil(current, dcu, turns, A, B, C)    
	%mi_addcircprop("circuitnumber", current,1=series|0=parallel)
	mi_addcircprop('Primary',       current,     1)
	mesh = dcu / 10;

	% draw windings
	mi_addnode(-B/2,0);
	mi_addnode(-B/2-A,0);
	mi_addnode(-B/2-A,C);
	mi_addnode(-B/2,C);
	mi_addnode(B/2,C);
	mi_addnode(B/2+A,0);
	mi_addnode(B/2+A,C);
	mi_addnode(B/2,0);

	mi_addsegment(-B/2,0,-B/2-A,0);
	mi_addsegment(-B/2-A,0,-B/2-A,C);
	mi_addsegment(-B/2-A,C,-B/2,C);
	mi_addsegment(-B/2,C,-B/2,0);
	mi_addsegment(B/2,0,B/2+A,0);
	mi_addsegment(B/2+A,0,B/2+A,C);
	mi_addsegment(B/2+A,C,B/2,C);
	mi_addsegment(B/2,C,B/2,0);

	%-- Copper Properties
	%-- PLUS
	mi_addblocklabel(B/2+A/2,C/3);
	mi_selectlabel(B/2+A/2,C/3);
	mi_setblockprop(sprintf('%0.1fmm',dcu),    1,       5 ,     'Primary', 0, 101, turns);
	mi_clearselected();
	%-- MINUS
	mi_addblocklabel(-B/2-A/2,C/3);
	mi_selectlabel(-B/2-A/2,C/3);
	mi_setblockprop(sprintf('%0.1fmm',dcu),    1,       5 ,     'Primary', 0, 101, -turns);
	mi_clearselected();
endfunction

function analyze(B, C)
	mi_saveas('.\\temp\\plots\\temp2.fem');
	mi_analyze();
	mi_loadsolution();
	mo_zoomnatural();
	mo_showdensityplot(-1,0,2.5,0,'mag')
	%mo_makeplot(PlotType,NumPoints,Filename,FileFormat)
	mo_addcontour(-B/2,C/2);
	mo_addcontour(B/2,C/2);
	mo_makeplot(1,300);
	mo_makeplot(1,300,'.\\temp\\numeric.txt',0)
	mi_clearselected();
endfunction
    
function calcInductance(current)    
	mo_groupselectblock(101);
	AJ1 = mo_blockintegral(0);
	BlockVolume = (mo_blockintegral(10));
	Pvx = mo_blockintegral(6);
	mo_clearblock();
	L11x = real(AJ1 / current / current);
	R11x = real(Pvx  / current / current) ;
	LL(1) = L11x;
	RR(1) = R11x;
	global inductance;
	global resistance;
	inductance = L11x*1e6/1000;
	resistance = R11x*1000;
endfunction

function printResults(inductance, resistance, current, turns, strands, dcu, A, B, C, Depth, Wcu)
	out = sprintf('Inductance..........: %0.3f mH',inductance);   disp(out);
	out = sprintf('Resistance..........: %0.0f mOhm',resistance);   disp(out);
	out = sprintf('Current.............: %0.0f A',current);   disp(out);
	out = sprintf('No of turns = %0.0f',turns);   disp(out);
	out = sprintf('%d x %0.2f mm',strands,dcu);   disp(out);
	out = sprintf('Width = %0.0f mm',2*A+B);   disp(out);
	out = sprintf('Height = %0.0f mm',C);   disp(out);
	out = sprintf('Depth = %0.0f mm',Depth);   disp(out);
	out = sprintf('Weight  = %0.2f kg',Wcu);   disp(out);
endfunction



initFEMM(Depth)
addMaterials(Cduct, strands, dcu)
setupWorkspace(A, B, C) 
drawCoil(current, dcu, turns, A, B, C)
analyze(B, C)
calcInductance(current)    
printResults(inductance, resistance, current, turns, strands, dcu, A, B, C, Depth, Wcu)

closefemm()

save ./temp/LandDCR.txt inductance resistance;
