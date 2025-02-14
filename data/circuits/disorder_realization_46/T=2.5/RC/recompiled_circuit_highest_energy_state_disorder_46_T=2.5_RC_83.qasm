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
rz(0.47705874) q[0];
sx q[0];
rz(2.8395489) q[0];
sx q[0];
rz(10.680847) q[0];
rz(0.37580252) q[1];
sx q[1];
rz(1.8355651) q[1];
sx q[1];
rz(10.553283) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9039584) q[0];
sx q[0];
rz(-2.4914143) q[0];
sx q[0];
rz(-2.0439434) q[0];
rz(-1.2641505) q[2];
sx q[2];
rz(-2.1696343) q[2];
sx q[2];
rz(1.5360338) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0998209) q[1];
sx q[1];
rz(-1.9536293) q[1];
sx q[1];
rz(1.0284996) q[1];
x q[2];
rz(-0.56803583) q[3];
sx q[3];
rz(-0.71692335) q[3];
sx q[3];
rz(0.64729819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6065373) q[2];
sx q[2];
rz(-0.29416072) q[2];
sx q[2];
rz(1.9425707) q[2];
rz(-2.8626275) q[3];
sx q[3];
rz(-1.6821034) q[3];
sx q[3];
rz(-3.103106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26175427) q[0];
sx q[0];
rz(-0.45670515) q[0];
sx q[0];
rz(-0.36619151) q[0];
rz(1.2884033) q[1];
sx q[1];
rz(-0.90694702) q[1];
sx q[1];
rz(-0.2624951) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85030445) q[0];
sx q[0];
rz(-1.8666963) q[0];
sx q[0];
rz(-1.759191) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3448703) q[2];
sx q[2];
rz(-2.886555) q[2];
sx q[2];
rz(-1.13305) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4387063) q[1];
sx q[1];
rz(-1.2186495) q[1];
sx q[1];
rz(-1.8789324) q[1];
rz(1.6672584) q[3];
sx q[3];
rz(-2.8854807) q[3];
sx q[3];
rz(-2.2369568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75258201) q[2];
sx q[2];
rz(-1.8730619) q[2];
sx q[2];
rz(2.9627723) q[2];
rz(3.1255152) q[3];
sx q[3];
rz(-2.6060846) q[3];
sx q[3];
rz(-2.2107562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43612424) q[0];
sx q[0];
rz(-0.12691623) q[0];
sx q[0];
rz(-2.5674168) q[0];
rz(0.15159675) q[1];
sx q[1];
rz(-0.44436374) q[1];
sx q[1];
rz(1.9134329) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3686217) q[0];
sx q[0];
rz(-0.54062343) q[0];
sx q[0];
rz(1.2631046) q[0];
x q[1];
rz(-0.068563383) q[2];
sx q[2];
rz(-2.1625963) q[2];
sx q[2];
rz(0.14321274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91635242) q[1];
sx q[1];
rz(-2.181899) q[1];
sx q[1];
rz(-3.0322187) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28801961) q[3];
sx q[3];
rz(-2.0374808) q[3];
sx q[3];
rz(-2.6096901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.14153081) q[2];
sx q[2];
rz(-2.5059293) q[2];
sx q[2];
rz(0.89111152) q[2];
rz(0.38854232) q[3];
sx q[3];
rz(-1.6308234) q[3];
sx q[3];
rz(2.8019606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163088) q[0];
sx q[0];
rz(-0.044965222) q[0];
sx q[0];
rz(-0.48211023) q[0];
rz(-0.66857839) q[1];
sx q[1];
rz(-1.6079638) q[1];
sx q[1];
rz(-0.63749981) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9264049) q[0];
sx q[0];
rz(-2.5026461) q[0];
sx q[0];
rz(-1.1730582) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89738557) q[2];
sx q[2];
rz(-1.4801133) q[2];
sx q[2];
rz(0.53340837) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.54957057) q[1];
sx q[1];
rz(-1.6046822) q[1];
sx q[1];
rz(3.0646044) q[1];
rz(2.8907475) q[3];
sx q[3];
rz(-1.9168087) q[3];
sx q[3];
rz(-3.0437247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1981226) q[2];
sx q[2];
rz(-2.1169457) q[2];
sx q[2];
rz(-0.060297273) q[2];
rz(0.30334011) q[3];
sx q[3];
rz(-3.0637432) q[3];
sx q[3];
rz(2.8764603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9918793) q[0];
sx q[0];
rz(-2.7880221) q[0];
sx q[0];
rz(2.7283332) q[0];
rz(-0.39516559) q[1];
sx q[1];
rz(-1.7839849) q[1];
sx q[1];
rz(-1.0023592) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7189594) q[0];
sx q[0];
rz(-2.0074685) q[0];
sx q[0];
rz(0.60499259) q[0];
x q[1];
rz(-3.1055519) q[2];
sx q[2];
rz(-1.3610351) q[2];
sx q[2];
rz(-2.710611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9262979) q[1];
sx q[1];
rz(-1.3560672) q[1];
sx q[1];
rz(-0.88820998) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2484364) q[3];
sx q[3];
rz(-1.3731706) q[3];
sx q[3];
rz(1.2569497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5137382) q[2];
sx q[2];
rz(-2.0444874) q[2];
sx q[2];
rz(1.6045137) q[2];
rz(1.9419935) q[3];
sx q[3];
rz(-2.0204085) q[3];
sx q[3];
rz(0.77170038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5143249) q[0];
sx q[0];
rz(-1.6384614) q[0];
sx q[0];
rz(2.7728873) q[0];
rz(0.74327028) q[1];
sx q[1];
rz(-2.5786046) q[1];
sx q[1];
rz(-0.35400131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1339851) q[0];
sx q[0];
rz(-2.2107894) q[0];
sx q[0];
rz(-0.065153393) q[0];
x q[1];
rz(2.6145489) q[2];
sx q[2];
rz(-1.8436925) q[2];
sx q[2];
rz(2.5986586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1917369) q[1];
sx q[1];
rz(-2.5402148) q[1];
sx q[1];
rz(-1.2527501) q[1];
rz(-pi) q[2];
rz(2.5052222) q[3];
sx q[3];
rz(-1.2546347) q[3];
sx q[3];
rz(3.0077601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1165498) q[2];
sx q[2];
rz(-1.9321238) q[2];
sx q[2];
rz(-0.0054736007) q[2];
rz(2.6239851) q[3];
sx q[3];
rz(-2.7766683) q[3];
sx q[3];
rz(-3.1143809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40590498) q[0];
sx q[0];
rz(-0.55957496) q[0];
sx q[0];
rz(2.9285808) q[0];
rz(-1.4951911) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(0.84943736) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6005858) q[0];
sx q[0];
rz(-1.6141506) q[0];
sx q[0];
rz(-0.61886532) q[0];
x q[1];
rz(-3.1108625) q[2];
sx q[2];
rz(-1.6140013) q[2];
sx q[2];
rz(1.2655366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35104677) q[1];
sx q[1];
rz(-2.4645269) q[1];
sx q[1];
rz(0.76404567) q[1];
rz(-pi) q[2];
rz(-0.15504239) q[3];
sx q[3];
rz(-2.0938805) q[3];
sx q[3];
rz(1.6301159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9269632) q[2];
sx q[2];
rz(-1.9081076) q[2];
sx q[2];
rz(-2.5978973) q[2];
rz(-0.26116535) q[3];
sx q[3];
rz(-1.5800913) q[3];
sx q[3];
rz(-0.15682596) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1866622) q[0];
sx q[0];
rz(-1.0533227) q[0];
sx q[0];
rz(2.9539811) q[0];
rz(-1.1232417) q[1];
sx q[1];
rz(-0.76485991) q[1];
sx q[1];
rz(0.16256464) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150972) q[0];
sx q[0];
rz(-2.2216068) q[0];
sx q[0];
rz(-1.9810173) q[0];
rz(-pi) q[1];
rz(-1.0467023) q[2];
sx q[2];
rz(-1.7225637) q[2];
sx q[2];
rz(-0.87226112) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9342207) q[1];
sx q[1];
rz(-1.5255063) q[1];
sx q[1];
rz(2.1846733) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1127693) q[3];
sx q[3];
rz(-0.9774607) q[3];
sx q[3];
rz(1.6360374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8673914) q[2];
sx q[2];
rz(-0.31495366) q[2];
sx q[2];
rz(-2.9401722) q[2];
rz(-2.6650186) q[3];
sx q[3];
rz(-1.3786992) q[3];
sx q[3];
rz(0.99982888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1808712) q[0];
sx q[0];
rz(-2.9044386) q[0];
sx q[0];
rz(2.5647105) q[0];
rz(0.30413973) q[1];
sx q[1];
rz(-2.0174593) q[1];
sx q[1];
rz(2.0374128) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7542412) q[0];
sx q[0];
rz(-1.1925444) q[0];
sx q[0];
rz(-3.0463112) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9530802) q[2];
sx q[2];
rz(-2.5645263) q[2];
sx q[2];
rz(0.12332502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6526281) q[1];
sx q[1];
rz(-1.9581989) q[1];
sx q[1];
rz(-2.1968958) q[1];
x q[2];
rz(3.125413) q[3];
sx q[3];
rz(-1.9619298) q[3];
sx q[3];
rz(-1.1333381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96678174) q[2];
sx q[2];
rz(-0.2468144) q[2];
sx q[2];
rz(-2.4583869) q[2];
rz(-1.5934058) q[3];
sx q[3];
rz(-2.4631409) q[3];
sx q[3];
rz(2.9233176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19112793) q[0];
sx q[0];
rz(-2.2948965) q[0];
sx q[0];
rz(-0.51093131) q[0];
rz(1.9966985) q[1];
sx q[1];
rz(-1.9547918) q[1];
sx q[1];
rz(1.5085545) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.78054) q[0];
sx q[0];
rz(-1.0298488) q[0];
sx q[0];
rz(1.279328) q[0];
x q[1];
rz(-1.0867003) q[2];
sx q[2];
rz(-1.3333313) q[2];
sx q[2];
rz(0.1456085) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3737693) q[1];
sx q[1];
rz(-1.1255742) q[1];
sx q[1];
rz(-1.7675464) q[1];
rz(-0.69798098) q[3];
sx q[3];
rz(-2.7578286) q[3];
sx q[3];
rz(2.021054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42023429) q[2];
sx q[2];
rz(-0.63174641) q[2];
sx q[2];
rz(-1.7101804) q[2];
rz(3.1259649) q[3];
sx q[3];
rz(-1.8764007) q[3];
sx q[3];
rz(2.6354852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2331727) q[0];
sx q[0];
rz(-1.40042) q[0];
sx q[0];
rz(-0.9216876) q[0];
rz(1.5064916) q[1];
sx q[1];
rz(-1.2163305) q[1];
sx q[1];
rz(-0.46179927) q[1];
rz(1.9533515) q[2];
sx q[2];
rz(-1.8391704) q[2];
sx q[2];
rz(3.0808629) q[2];
rz(-1.9962068) q[3];
sx q[3];
rz(-0.085193188) q[3];
sx q[3];
rz(-1.9151158) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
