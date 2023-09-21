OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(-1.707466) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(-2.5399962) q[1];
sx q[1];
rz(-0.4184202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7941064) q[0];
sx q[0];
rz(-2.5972164) q[0];
sx q[0];
rz(2.0967336) q[0];
rz(1.4346052) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(1.1510804) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5601555) q[1];
sx q[1];
rz(-1.1168224) q[1];
sx q[1];
rz(-0.27172383) q[1];
rz(-2.0099785) q[3];
sx q[3];
rz(-1.3703128) q[3];
sx q[3];
rz(0.85103121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.819954) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-0.83797541) q[2];
rz(2.6485802) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74719602) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.2778506) q[0];
rz(0.17678075) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(0.4321672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5837001) q[0];
sx q[0];
rz(-2.0895134) q[0];
sx q[0];
rz(-1.2516663) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.400153) q[2];
sx q[2];
rz(-1.1740985) q[2];
sx q[2];
rz(0.19043365) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1951616) q[1];
sx q[1];
rz(-2.6841607) q[1];
sx q[1];
rz(1.3034348) q[1];
rz(0.024859419) q[3];
sx q[3];
rz(-1.0059352) q[3];
sx q[3];
rz(0.75627518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5923578) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(-2.6548927) q[2];
rz(1.3782079) q[3];
sx q[3];
rz(-1.8816201) q[3];
sx q[3];
rz(-2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(2.7242463) q[0];
rz(1.4886645) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(-2.6352077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2877809) q[0];
sx q[0];
rz(-1.3940485) q[0];
sx q[0];
rz(1.9275097) q[0];
rz(-pi) q[1];
rz(-0.73451368) q[2];
sx q[2];
rz(-1.3909512) q[2];
sx q[2];
rz(-2.432446) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.312078) q[1];
sx q[1];
rz(-1.2631589) q[1];
sx q[1];
rz(-2.679146) q[1];
x q[2];
rz(-1.1261602) q[3];
sx q[3];
rz(-0.4053084) q[3];
sx q[3];
rz(-2.0733548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4425519) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(-0.59147269) q[2];
rz(-2.5555723) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(-1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9716924) q[0];
sx q[0];
rz(-0.62830347) q[0];
sx q[0];
rz(-2.3024094) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(1.5930088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65028134) q[0];
sx q[0];
rz(-2.6822753) q[0];
sx q[0];
rz(-1.7772872) q[0];
rz(-pi) q[1];
rz(-2.2976539) q[2];
sx q[2];
rz(-2.2288449) q[2];
sx q[2];
rz(2.7968) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81739391) q[1];
sx q[1];
rz(-0.84016582) q[1];
sx q[1];
rz(-0.52449951) q[1];
x q[2];
rz(1.8539092) q[3];
sx q[3];
rz(-1.9861756) q[3];
sx q[3];
rz(2.0302041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30535355) q[2];
sx q[2];
rz(-2.2576015) q[2];
sx q[2];
rz(-3.0419066) q[2];
rz(0.95885197) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(-1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754159) q[0];
sx q[0];
rz(-1.7594936) q[0];
sx q[0];
rz(-0.25594041) q[0];
rz(0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(2.3815313) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3292023) q[0];
sx q[0];
rz(-2.985552) q[0];
sx q[0];
rz(0.72326707) q[0];
rz(1.4270093) q[2];
sx q[2];
rz(-2.723105) q[2];
sx q[2];
rz(-0.33868044) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8885986) q[1];
sx q[1];
rz(-0.18054403) q[1];
sx q[1];
rz(-0.79043364) q[1];
x q[2];
rz(2.2171668) q[3];
sx q[3];
rz(-1.8674208) q[3];
sx q[3];
rz(-2.7588206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(0.6742397) q[2];
rz(-0.21480602) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(-0.017344346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53133416) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(-1.0669605) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63283352) q[0];
sx q[0];
rz(-1.5914704) q[0];
sx q[0];
rz(1.585292) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42491575) q[2];
sx q[2];
rz(-1.1669461) q[2];
sx q[2];
rz(-2.0770819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8384335) q[1];
sx q[1];
rz(-1.618297) q[1];
sx q[1];
rz(0.82364239) q[1];
rz(0.79640572) q[3];
sx q[3];
rz(-1.685433) q[3];
sx q[3];
rz(0.20533268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.431488) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(2.5816494) q[2];
rz(-2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(-0.30549756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840435) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(-1.7204826) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(-2.2834159) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2150543) q[0];
sx q[0];
rz(-0.19556043) q[0];
sx q[0];
rz(2.2901448) q[0];
x q[1];
rz(-3.0965205) q[2];
sx q[2];
rz(-1.1606693) q[2];
sx q[2];
rz(-0.28442597) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9550025) q[1];
sx q[1];
rz(-1.5555256) q[1];
sx q[1];
rz(-3.112622) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51123294) q[3];
sx q[3];
rz(-0.63496642) q[3];
sx q[3];
rz(0.2583897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28785607) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(1.8161592) q[2];
rz(-0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(-0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778075) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(2.7668787) q[0];
rz(-0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6357248) q[0];
sx q[0];
rz(-0.74698193) q[0];
sx q[0];
rz(2.5542459) q[0];
rz(-pi) q[1];
rz(1.5937514) q[2];
sx q[2];
rz(-2.288753) q[2];
sx q[2];
rz(0.1644451) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0950349) q[1];
sx q[1];
rz(-1.7525502) q[1];
sx q[1];
rz(-0.05038105) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3091062) q[3];
sx q[3];
rz(-1.040254) q[3];
sx q[3];
rz(1.4025276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4902041) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(-0.46869579) q[2];
rz(1.9474585) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(1.2619031) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(1.6040241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3106829) q[0];
sx q[0];
rz(-1.791782) q[0];
sx q[0];
rz(1.2587147) q[0];
x q[1];
rz(1.4666918) q[2];
sx q[2];
rz(-1.1541919) q[2];
sx q[2];
rz(-0.57550752) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9201587) q[1];
sx q[1];
rz(-2.4281574) q[1];
sx q[1];
rz(-0.17381298) q[1];
rz(-pi) q[2];
rz(-0.92026199) q[3];
sx q[3];
rz(-2.5839621) q[3];
sx q[3];
rz(0.44616163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1016772) q[2];
sx q[2];
rz(-0.95255178) q[2];
sx q[2];
rz(-0.76134479) q[2];
rz(0.90041655) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(0.049023978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(-2.9246869) q[0];
rz(0.63198173) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(-0.95473081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023708658) q[0];
sx q[0];
rz(-2.3291991) q[0];
sx q[0];
rz(-0.448416) q[0];
rz(2.5920792) q[2];
sx q[2];
rz(-1.9351442) q[2];
sx q[2];
rz(3.0669616) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.875647) q[1];
sx q[1];
rz(-1.4946283) q[1];
sx q[1];
rz(-0.45653685) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93654376) q[3];
sx q[3];
rz(-1.9359971) q[3];
sx q[3];
rz(1.6983502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1251936) q[2];
sx q[2];
rz(-1.23896) q[2];
sx q[2];
rz(-0.94474244) q[2];
rz(-2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(-2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5007297) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(0.62467081) q[2];
sx q[2];
rz(-2.2041337) q[2];
sx q[2];
rz(0.49908257) q[2];
rz(-1.1273884) q[3];
sx q[3];
rz(-1.3848806) q[3];
sx q[3];
rz(1.0564907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
