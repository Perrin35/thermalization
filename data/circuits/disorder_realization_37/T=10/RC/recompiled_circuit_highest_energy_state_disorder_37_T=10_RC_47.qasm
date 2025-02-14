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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(-2.6824644) q[0];
rz(-0.65161172) q[1];
sx q[1];
rz(-1.7348644) q[1];
sx q[1];
rz(-2.9098517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15329862) q[0];
sx q[0];
rz(-2.7911515) q[0];
sx q[0];
rz(-2.9399583) q[0];
rz(0.54643537) q[2];
sx q[2];
rz(-1.6197512) q[2];
sx q[2];
rz(-0.10412439) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5197088) q[1];
sx q[1];
rz(-2.1823898) q[1];
sx q[1];
rz(-1.1033464) q[1];
x q[2];
rz(1.5292589) q[3];
sx q[3];
rz(-1.9088863) q[3];
sx q[3];
rz(1.5134144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.090652466) q[2];
sx q[2];
rz(-0.53477627) q[2];
sx q[2];
rz(-1.862662) q[2];
rz(-2.2517962) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(2.8029627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5559674) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(-2.2157748) q[0];
rz(-1.0649118) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(3.056114) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8166072) q[0];
sx q[0];
rz(-2.7517635) q[0];
sx q[0];
rz(3.1329324) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34833724) q[2];
sx q[2];
rz(-2.1536963) q[2];
sx q[2];
rz(-0.23886853) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6798415) q[1];
sx q[1];
rz(-1.0437168) q[1];
sx q[1];
rz(2.8746469) q[1];
x q[2];
rz(1.1033589) q[3];
sx q[3];
rz(-0.93884838) q[3];
sx q[3];
rz(1.8935209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.815879) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(-2.7549287) q[2];
rz(1.9246842) q[3];
sx q[3];
rz(-0.34438008) q[3];
sx q[3];
rz(-0.2200505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81095186) q[0];
sx q[0];
rz(-2.7702259) q[0];
sx q[0];
rz(-0.49705848) q[0];
rz(-2.0992384) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(-0.29464468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9884856) q[0];
sx q[0];
rz(-1.238958) q[0];
sx q[0];
rz(1.6436623) q[0];
x q[1];
rz(-2.5025236) q[2];
sx q[2];
rz(-1.7099755) q[2];
sx q[2];
rz(-0.39060171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4528447) q[1];
sx q[1];
rz(-2.0032681) q[1];
sx q[1];
rz(-0.8853064) q[1];
rz(-pi) q[2];
rz(-0.96986846) q[3];
sx q[3];
rz(-1.4517541) q[3];
sx q[3];
rz(2.5210019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3543388) q[2];
sx q[2];
rz(-0.86091176) q[2];
sx q[2];
rz(1.5617237) q[2];
rz(1.076237) q[3];
sx q[3];
rz(-1.302224) q[3];
sx q[3];
rz(-0.92897433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314986) q[0];
sx q[0];
rz(-1.6153233) q[0];
sx q[0];
rz(-0.22931799) q[0];
rz(1.5062821) q[1];
sx q[1];
rz(-1.458026) q[1];
sx q[1];
rz(-3.07952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784883) q[0];
sx q[0];
rz(-2.0651428) q[0];
sx q[0];
rz(0.93017471) q[0];
rz(-1.1787291) q[2];
sx q[2];
rz(-1.829756) q[2];
sx q[2];
rz(-0.94253892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2885372) q[1];
sx q[1];
rz(-1.2934577) q[1];
sx q[1];
rz(0.71373516) q[1];
rz(-pi) q[2];
rz(-1.6639978) q[3];
sx q[3];
rz(-1.3706932) q[3];
sx q[3];
rz(-2.4148108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.092992358) q[2];
sx q[2];
rz(-2.2255662) q[2];
sx q[2];
rz(-1.830706) q[2];
rz(-3.1033031) q[3];
sx q[3];
rz(-1.7891276) q[3];
sx q[3];
rz(-0.37461764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49486092) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(0.89299655) q[0];
rz(2.112174) q[1];
sx q[1];
rz(-2.4118377) q[1];
sx q[1];
rz(3.0457048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9150118) q[0];
sx q[0];
rz(-1.9862729) q[0];
sx q[0];
rz(1.0092495) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4818727) q[2];
sx q[2];
rz(-1.6409573) q[2];
sx q[2];
rz(-1.905575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.694987) q[1];
sx q[1];
rz(-2.1220653) q[1];
sx q[1];
rz(2.0010082) q[1];
rz(-pi) q[2];
rz(2.5378102) q[3];
sx q[3];
rz(-2.5386435) q[3];
sx q[3];
rz(-0.83151885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9868077) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(0.42243877) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(-2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4130037) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(3.0244306) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(-2.07043) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9566222) q[0];
sx q[0];
rz(-1.8529467) q[0];
sx q[0];
rz(2.6376702) q[0];
x q[1];
rz(-0.44354673) q[2];
sx q[2];
rz(-2.5534592) q[2];
sx q[2];
rz(0.085970446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92381645) q[1];
sx q[1];
rz(-0.72715302) q[1];
sx q[1];
rz(-1.9132861) q[1];
rz(0.31757327) q[3];
sx q[3];
rz(-1.4151207) q[3];
sx q[3];
rz(-2.6359216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.57227197) q[2];
sx q[2];
rz(-2.5666777) q[2];
sx q[2];
rz(3.102109) q[2];
rz(-3.0651921) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(2.6881645) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4385248) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(-1.9901265) q[1];
sx q[1];
rz(-1.6116424) q[1];
sx q[1];
rz(-1.8340402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3261953) q[0];
sx q[0];
rz(-1.2457677) q[0];
sx q[0];
rz(-2.9103878) q[0];
rz(0.18270638) q[2];
sx q[2];
rz(-2.7596843) q[2];
sx q[2];
rz(-2.4289102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3486191) q[1];
sx q[1];
rz(-1.359531) q[1];
sx q[1];
rz(-2.0729077) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3780955) q[3];
sx q[3];
rz(-2.5268838) q[3];
sx q[3];
rz(-2.4845882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82039708) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(-2.5879228) q[2];
rz(0.6959483) q[3];
sx q[3];
rz(-1.4724933) q[3];
sx q[3];
rz(-1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11947908) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(-2.8346862) q[0];
rz(-1.8652929) q[1];
sx q[1];
rz(-0.46638322) q[1];
sx q[1];
rz(-1.2409522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0690445) q[0];
sx q[0];
rz(-0.47397754) q[0];
sx q[0];
rz(1.2233673) q[0];
x q[1];
rz(2.5436756) q[2];
sx q[2];
rz(-1.2094524) q[2];
sx q[2];
rz(1.8449699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3731662) q[1];
sx q[1];
rz(-1.2618999) q[1];
sx q[1];
rz(-1.373148) q[1];
rz(-pi) q[2];
rz(-2.5316174) q[3];
sx q[3];
rz(-1.9053359) q[3];
sx q[3];
rz(1.1776678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3500195) q[2];
sx q[2];
rz(-0.80524033) q[2];
sx q[2];
rz(1.2903068) q[2];
rz(-0.6811412) q[3];
sx q[3];
rz(-0.9001503) q[3];
sx q[3];
rz(-1.5017989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8660368) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(-1.9679029) q[0];
rz(2.5545919) q[1];
sx q[1];
rz(-2.360011) q[1];
sx q[1];
rz(0.079708286) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82725924) q[0];
sx q[0];
rz(-1.9406576) q[0];
sx q[0];
rz(-2.8989559) q[0];
rz(-pi) q[1];
rz(-0.22144145) q[2];
sx q[2];
rz(-1.4613073) q[2];
sx q[2];
rz(1.2564645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6708128) q[1];
sx q[1];
rz(-2.1683886) q[1];
sx q[1];
rz(-2.2714991) q[1];
rz(1.2740785) q[3];
sx q[3];
rz(-2.35459) q[3];
sx q[3];
rz(-2.3001461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6012663) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(1.8076753) q[2];
rz(-1.5506844) q[3];
sx q[3];
rz(-1.0100789) q[3];
sx q[3];
rz(-0.18868119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48225668) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(-1.851086) q[0];
rz(-0.036529649) q[1];
sx q[1];
rz(-1.1643658) q[1];
sx q[1];
rz(-1.0443002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8528906) q[0];
sx q[0];
rz(-1.9963309) q[0];
sx q[0];
rz(0.3279) q[0];
rz(-pi) q[1];
rz(0.26728435) q[2];
sx q[2];
rz(-1.2517831) q[2];
sx q[2];
rz(0.50179447) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3022233) q[1];
sx q[1];
rz(-1.7249134) q[1];
sx q[1];
rz(-0.83854143) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4560111) q[3];
sx q[3];
rz(-1.8978531) q[3];
sx q[3];
rz(2.9714874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9145987) q[2];
sx q[2];
rz(-1.031216) q[2];
sx q[2];
rz(-1.3247789) q[2];
rz(-0.75657183) q[3];
sx q[3];
rz(-2.2533267) q[3];
sx q[3];
rz(-2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4751547) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(1.4229763) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(-0.46427609) q[2];
sx q[2];
rz(-0.79002476) q[2];
sx q[2];
rz(2.644047) q[2];
rz(2.4800469) q[3];
sx q[3];
rz(-1.9521803) q[3];
sx q[3];
rz(3.1357756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
