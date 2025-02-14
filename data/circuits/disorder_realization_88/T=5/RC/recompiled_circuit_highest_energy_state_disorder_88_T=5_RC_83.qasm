OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54230827) q[0];
sx q[0];
rz(-0.13442726) q[0];
sx q[0];
rz(2.0943213) q[0];
rz(-0.32416999) q[1];
sx q[1];
rz(-0.14961641) q[1];
sx q[1];
rz(-2.1967998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1143998) q[0];
sx q[0];
rz(-1.8437705) q[0];
sx q[0];
rz(1.9045715) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9247577) q[2];
sx q[2];
rz(-1.3831209) q[2];
sx q[2];
rz(-3.0800703) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.289669) q[1];
sx q[1];
rz(-2.0229368) q[1];
sx q[1];
rz(0.55638066) q[1];
rz(1.0881926) q[3];
sx q[3];
rz(-1.6783829) q[3];
sx q[3];
rz(1.7501329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85113877) q[2];
sx q[2];
rz(-1.4522499) q[2];
sx q[2];
rz(1.5770844) q[2];
rz(0.56646281) q[3];
sx q[3];
rz(-0.62140673) q[3];
sx q[3];
rz(-1.8433146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31747776) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(-2.4216006) q[0];
rz(0.27436817) q[1];
sx q[1];
rz(-0.73148483) q[1];
sx q[1];
rz(2.5879587) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047129) q[0];
sx q[0];
rz(-1.9229445) q[0];
sx q[0];
rz(3.0180305) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4963412) q[2];
sx q[2];
rz(-2.650839) q[2];
sx q[2];
rz(3.0440999) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.28277507) q[1];
sx q[1];
rz(-1.2056279) q[1];
sx q[1];
rz(-2.3479985) q[1];
x q[2];
rz(0.7577866) q[3];
sx q[3];
rz(-0.67000853) q[3];
sx q[3];
rz(2.3520833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9686034) q[2];
sx q[2];
rz(-1.9393238) q[2];
sx q[2];
rz(2.7776862) q[2];
rz(-3.0299752) q[3];
sx q[3];
rz(-0.93897396) q[3];
sx q[3];
rz(-2.3811471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97742057) q[0];
sx q[0];
rz(-2.3168679) q[0];
sx q[0];
rz(0.0053996276) q[0];
rz(-2.3258356) q[1];
sx q[1];
rz(-2.0033671) q[1];
sx q[1];
rz(2.7556748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513228) q[0];
sx q[0];
rz(-1.9093336) q[0];
sx q[0];
rz(1.2744034) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0089985) q[2];
sx q[2];
rz(-1.2616065) q[2];
sx q[2];
rz(-2.5872292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93497998) q[1];
sx q[1];
rz(-0.40217186) q[1];
sx q[1];
rz(1.9002756) q[1];
rz(-pi) q[2];
rz(0.13428646) q[3];
sx q[3];
rz(-1.9818174) q[3];
sx q[3];
rz(-2.3255796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58891121) q[2];
sx q[2];
rz(-0.83325714) q[2];
sx q[2];
rz(1.3820648) q[2];
rz(2.9803993) q[3];
sx q[3];
rz(-1.4313982) q[3];
sx q[3];
rz(2.5126357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2564119) q[0];
sx q[0];
rz(-1.2762524) q[0];
sx q[0];
rz(-1.3389583) q[0];
rz(1.3171875) q[1];
sx q[1];
rz(-0.97511292) q[1];
sx q[1];
rz(-0.29979527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6000344) q[0];
sx q[0];
rz(-2.2081714) q[0];
sx q[0];
rz(-1.1336898) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6671403) q[2];
sx q[2];
rz(-1.0336733) q[2];
sx q[2];
rz(-0.65309292) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2721709) q[1];
sx q[1];
rz(-1.7079087) q[1];
sx q[1];
rz(1.2161944) q[1];
x q[2];
rz(2.2065877) q[3];
sx q[3];
rz(-1.6775369) q[3];
sx q[3];
rz(-0.028117953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51912159) q[2];
sx q[2];
rz(-2.2032479) q[2];
sx q[2];
rz(0.31912121) q[2];
rz(0.090911344) q[3];
sx q[3];
rz(-0.14486434) q[3];
sx q[3];
rz(1.7849785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262481) q[0];
sx q[0];
rz(-2.321796) q[0];
sx q[0];
rz(3.1392642) q[0];
rz(-1.5709411) q[1];
sx q[1];
rz(-2.3410773) q[1];
sx q[1];
rz(-1.194582) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9228421) q[0];
sx q[0];
rz(-0.19521579) q[0];
sx q[0];
rz(0.65252916) q[0];
rz(-3.1011337) q[2];
sx q[2];
rz(-0.85340696) q[2];
sx q[2];
rz(2.3641912) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4264226) q[1];
sx q[1];
rz(-1.914901) q[1];
sx q[1];
rz(0.10870966) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8525328) q[3];
sx q[3];
rz(-1.0322555) q[3];
sx q[3];
rz(2.0592214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.883256) q[2];
sx q[2];
rz(-0.77815762) q[2];
sx q[2];
rz(0.077433132) q[2];
rz(2.1954913) q[3];
sx q[3];
rz(-0.48282048) q[3];
sx q[3];
rz(-2.6879123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2133863) q[0];
sx q[0];
rz(-0.086406924) q[0];
sx q[0];
rz(1.4148096) q[0];
rz(2.4746223) q[1];
sx q[1];
rz(-1.7844776) q[1];
sx q[1];
rz(2.736843) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0694097) q[0];
sx q[0];
rz(-1.8526005) q[0];
sx q[0];
rz(0.3148766) q[0];
x q[1];
rz(2.8799492) q[2];
sx q[2];
rz(-2.2115797) q[2];
sx q[2];
rz(1.6710168) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.09913832) q[1];
sx q[1];
rz(-2.6746866) q[1];
sx q[1];
rz(1.7068638) q[1];
rz(2.7795198) q[3];
sx q[3];
rz(-1.1112461) q[3];
sx q[3];
rz(-1.5611695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3177967) q[2];
sx q[2];
rz(-0.78936374) q[2];
sx q[2];
rz(-0.44343597) q[2];
rz(-0.41380841) q[3];
sx q[3];
rz(-0.2897073) q[3];
sx q[3];
rz(-2.2558291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39171788) q[0];
sx q[0];
rz(-1.3626008) q[0];
sx q[0];
rz(0.66746563) q[0];
rz(-0.04774566) q[1];
sx q[1];
rz(-2.0011438) q[1];
sx q[1];
rz(-1.1183636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3453433) q[0];
sx q[0];
rz(-2.44758) q[0];
sx q[0];
rz(-0.66008499) q[0];
rz(-3.0202853) q[2];
sx q[2];
rz(-2.1459142) q[2];
sx q[2];
rz(0.93851575) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27208313) q[1];
sx q[1];
rz(-2.1456302) q[1];
sx q[1];
rz(1.065567) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0712549) q[3];
sx q[3];
rz(-1.5111672) q[3];
sx q[3];
rz(-1.3524551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.095801) q[2];
sx q[2];
rz(-1.8818776) q[2];
sx q[2];
rz(0.0087139159) q[2];
rz(1.2039315) q[3];
sx q[3];
rz(-0.34398505) q[3];
sx q[3];
rz(0.022139942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1190204) q[0];
sx q[0];
rz(-2.75596) q[0];
sx q[0];
rz(-0.88985306) q[0];
rz(1.285137) q[1];
sx q[1];
rz(-2.0882008) q[1];
sx q[1];
rz(3.0056675) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21273577) q[0];
sx q[0];
rz(-1.341916) q[0];
sx q[0];
rz(1.5221734) q[0];
x q[1];
rz(-2.459002) q[2];
sx q[2];
rz(-0.75645489) q[2];
sx q[2];
rz(-0.83908778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4299851) q[1];
sx q[1];
rz(-1.2384336) q[1];
sx q[1];
rz(1.3620939) q[1];
x q[2];
rz(-0.14304026) q[3];
sx q[3];
rz(-1.2238127) q[3];
sx q[3];
rz(-8/(1*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4385684) q[2];
sx q[2];
rz(-0.88006222) q[2];
sx q[2];
rz(-1.0620037) q[2];
rz(-2.2117173) q[3];
sx q[3];
rz(-0.78571856) q[3];
sx q[3];
rz(2.3532531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36282614) q[0];
sx q[0];
rz(-2.7500948) q[0];
sx q[0];
rz(2.7407001) q[0];
rz(2.3800384) q[1];
sx q[1];
rz(-1.4925894) q[1];
sx q[1];
rz(-2.3525995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61643314) q[0];
sx q[0];
rz(-0.1285006) q[0];
sx q[0];
rz(-0.12501053) q[0];
x q[1];
rz(-2.2458399) q[2];
sx q[2];
rz(-1.6924685) q[2];
sx q[2];
rz(2.3020803) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6500864) q[1];
sx q[1];
rz(-1.6900926) q[1];
sx q[1];
rz(0.59081282) q[1];
rz(1.5212359) q[3];
sx q[3];
rz(-1.0465099) q[3];
sx q[3];
rz(1.4567284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.57364982) q[2];
sx q[2];
rz(-1.2929792) q[2];
sx q[2];
rz(0.72081494) q[2];
rz(1.6503664) q[3];
sx q[3];
rz(-0.78250116) q[3];
sx q[3];
rz(-1.8318374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1184621) q[0];
sx q[0];
rz(-0.11496249) q[0];
sx q[0];
rz(0.7777099) q[0];
rz(2.481781) q[1];
sx q[1];
rz(-0.86645627) q[1];
sx q[1];
rz(-2.7515817) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56905079) q[0];
sx q[0];
rz(-0.13992289) q[0];
sx q[0];
rz(-1.8403017) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71163746) q[2];
sx q[2];
rz(-2.2588428) q[2];
sx q[2];
rz(-1.2104863) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.240814) q[1];
sx q[1];
rz(-1.972318) q[1];
sx q[1];
rz(2.8177849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50647363) q[3];
sx q[3];
rz(-2.1586543) q[3];
sx q[3];
rz(1.7755058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68221349) q[2];
sx q[2];
rz(-2.7028658) q[2];
sx q[2];
rz(0.48975804) q[2];
rz(-2.8343936) q[3];
sx q[3];
rz(-0.92370737) q[3];
sx q[3];
rz(-1.7507621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3601892) q[0];
sx q[0];
rz(-1.4908981) q[0];
sx q[0];
rz(1.4733462) q[0];
rz(-2.0669943) q[1];
sx q[1];
rz(-2.2886724) q[1];
sx q[1];
rz(1.2235175) q[1];
rz(-2.3831153) q[2];
sx q[2];
rz(-0.52634326) q[2];
sx q[2];
rz(1.1545622) q[2];
rz(-2.1297602) q[3];
sx q[3];
rz(-1.0444416) q[3];
sx q[3];
rz(-2.039644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
