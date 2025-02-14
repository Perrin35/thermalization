OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45159856) q[0];
sx q[0];
rz(-0.30878433) q[0];
sx q[0];
rz(-0.2398332) q[0];
rz(2.580515) q[1];
sx q[1];
rz(-0.74220389) q[1];
sx q[1];
rz(-1.6286558) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0081019199) q[0];
sx q[0];
rz(-1.9372131) q[0];
sx q[0];
rz(-3.0045463) q[0];
x q[1];
rz(-2.0530897) q[2];
sx q[2];
rz(-2.6183207) q[2];
sx q[2];
rz(2.5918617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7907432) q[1];
sx q[1];
rz(-0.51209699) q[1];
sx q[1];
rz(-2.8382906) q[1];
rz(-pi) q[2];
rz(3.0835864) q[3];
sx q[3];
rz(-1.4027624) q[3];
sx q[3];
rz(1.9498273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1538887) q[2];
sx q[2];
rz(-2.14553) q[2];
sx q[2];
rz(-1.055701) q[2];
rz(-0.40134564) q[3];
sx q[3];
rz(-1.515712) q[3];
sx q[3];
rz(-1.5041941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6317247) q[0];
sx q[0];
rz(-2.8724176) q[0];
sx q[0];
rz(1.4058231) q[0];
rz(0.74613219) q[1];
sx q[1];
rz(-1.1765307) q[1];
sx q[1];
rz(-3.0768118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1593744) q[0];
sx q[0];
rz(-1.9663334) q[0];
sx q[0];
rz(-0.48898029) q[0];
rz(-2.9322036) q[2];
sx q[2];
rz(-2.4509666) q[2];
sx q[2];
rz(-2.2356981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3105615) q[1];
sx q[1];
rz(-1.4906724) q[1];
sx q[1];
rz(-1.4270272) q[1];
x q[2];
rz(1.359932) q[3];
sx q[3];
rz(-1.5393889) q[3];
sx q[3];
rz(-1.5109911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.259321) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(0.041291324) q[2];
rz(1.1193554) q[3];
sx q[3];
rz(-1.7740403) q[3];
sx q[3];
rz(-2.0004499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.816788) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(-0.69931716) q[0];
rz(-2.518867) q[1];
sx q[1];
rz(-0.8131578) q[1];
sx q[1];
rz(-0.25845382) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1656483) q[0];
sx q[0];
rz(-1.5585525) q[0];
sx q[0];
rz(-0.6338288) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7202366) q[2];
sx q[2];
rz(-2.4838313) q[2];
sx q[2];
rz(2.2542816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.099951) q[1];
sx q[1];
rz(-1.1219684) q[1];
sx q[1];
rz(-1.5931409) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4914037) q[3];
sx q[3];
rz(-0.43118048) q[3];
sx q[3];
rz(-0.10583793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2059325) q[2];
sx q[2];
rz(-1.6814597) q[2];
sx q[2];
rz(-0.70651954) q[2];
rz(-0.62890729) q[3];
sx q[3];
rz(-2.3157412) q[3];
sx q[3];
rz(2.0188873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0141456) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(-2.1674147) q[0];
rz(1.6575419) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(-1.0008224) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38353048) q[0];
sx q[0];
rz(-2.5208726) q[0];
sx q[0];
rz(-1.8412983) q[0];
rz(-pi) q[1];
rz(0.77299574) q[2];
sx q[2];
rz(-2.2064035) q[2];
sx q[2];
rz(2.5781812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.448424) q[1];
sx q[1];
rz(-2.1285372) q[1];
sx q[1];
rz(-1.0156469) q[1];
x q[2];
rz(-2.820904) q[3];
sx q[3];
rz(-1.6178035) q[3];
sx q[3];
rz(-0.21002029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68690825) q[2];
sx q[2];
rz(-1.800622) q[2];
sx q[2];
rz(-0.81721133) q[2];
rz(-0.59598437) q[3];
sx q[3];
rz(-1.417336) q[3];
sx q[3];
rz(1.3982841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3572094) q[0];
sx q[0];
rz(-1.7359808) q[0];
sx q[0];
rz(1.0719517) q[0];
rz(-2.0042073) q[1];
sx q[1];
rz(-1.6268566) q[1];
sx q[1];
rz(2.9683108) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6023077) q[0];
sx q[0];
rz(-2.0208911) q[0];
sx q[0];
rz(1.1316687) q[0];
x q[1];
rz(-2.7543147) q[2];
sx q[2];
rz(-2.0301295) q[2];
sx q[2];
rz(-1.1371374) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73663354) q[1];
sx q[1];
rz(-1.6425321) q[1];
sx q[1];
rz(2.328675) q[1];
rz(-2.9786879) q[3];
sx q[3];
rz(-1.6626193) q[3];
sx q[3];
rz(-3.0260928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5002354) q[2];
sx q[2];
rz(-2.6884029) q[2];
sx q[2];
rz(0.87453169) q[2];
rz(0.66679653) q[3];
sx q[3];
rz(-1.4614481) q[3];
sx q[3];
rz(0.82505208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2905529) q[0];
sx q[0];
rz(-0.14851004) q[0];
sx q[0];
rz(-2.9669951) q[0];
rz(-2.5945276) q[1];
sx q[1];
rz(-0.51536307) q[1];
sx q[1];
rz(-1.2778767) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1924612) q[0];
sx q[0];
rz(-1.2843848) q[0];
sx q[0];
rz(-0.064563036) q[0];
rz(2.6447884) q[2];
sx q[2];
rz(-0.81648177) q[2];
sx q[2];
rz(-1.4085438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7050347) q[1];
sx q[1];
rz(-0.83022699) q[1];
sx q[1];
rz(0.48506801) q[1];
rz(-pi) q[2];
rz(-2.296769) q[3];
sx q[3];
rz(-2.3298744) q[3];
sx q[3];
rz(-0.38540781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7586907) q[2];
sx q[2];
rz(-2.8768657) q[2];
sx q[2];
rz(0.36941377) q[2];
rz(2.3139125) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(0.15414342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5740042) q[0];
sx q[0];
rz(-0.91667691) q[0];
sx q[0];
rz(0.23707238) q[0];
rz(-1.392662) q[1];
sx q[1];
rz(-1.1940414) q[1];
sx q[1];
rz(-2.1377835) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7866681) q[0];
sx q[0];
rz(-2.0046429) q[0];
sx q[0];
rz(1.3870441) q[0];
rz(-pi) q[1];
rz(1.1091484) q[2];
sx q[2];
rz(-0.91141111) q[2];
sx q[2];
rz(-2.387799) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4752581) q[1];
sx q[1];
rz(-1.043519) q[1];
sx q[1];
rz(1.9557245) q[1];
rz(-2.5828894) q[3];
sx q[3];
rz(-1.1538343) q[3];
sx q[3];
rz(0.95181634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6134593) q[2];
sx q[2];
rz(-1.5626835) q[2];
sx q[2];
rz(-2.6252739) q[2];
rz(2.3033219) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(-2.9576438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8539921) q[0];
sx q[0];
rz(-2.8599399) q[0];
sx q[0];
rz(-0.91424346) q[0];
rz(0.5842579) q[1];
sx q[1];
rz(-1.7353568) q[1];
sx q[1];
rz(2.2343238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4582918) q[0];
sx q[0];
rz(-2.2193546) q[0];
sx q[0];
rz(-0.022819937) q[0];
x q[1];
rz(-1.3970988) q[2];
sx q[2];
rz(-0.96250421) q[2];
sx q[2];
rz(-2.4243958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71675352) q[1];
sx q[1];
rz(-2.3100634) q[1];
sx q[1];
rz(0.6296954) q[1];
rz(-pi) q[2];
rz(2.8738638) q[3];
sx q[3];
rz(-2.1526823) q[3];
sx q[3];
rz(-1.3022334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82627901) q[2];
sx q[2];
rz(-1.7118914) q[2];
sx q[2];
rz(0.86177525) q[2];
rz(1.2365384) q[3];
sx q[3];
rz(-3.0573513) q[3];
sx q[3];
rz(1.2110565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64060193) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(-1.030141) q[0];
rz(-2.8252699) q[1];
sx q[1];
rz(-1.7681237) q[1];
sx q[1];
rz(2.8533459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.41053) q[0];
sx q[0];
rz(-1.9707435) q[0];
sx q[0];
rz(1.7284786) q[0];
rz(-2.5000942) q[2];
sx q[2];
rz(-2.3833123) q[2];
sx q[2];
rz(-1.7816133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32375755) q[1];
sx q[1];
rz(-0.97684469) q[1];
sx q[1];
rz(-2.6282534) q[1];
x q[2];
rz(-1.3449617) q[3];
sx q[3];
rz(-0.91167456) q[3];
sx q[3];
rz(-2.3171605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2890702) q[2];
sx q[2];
rz(-1.728629) q[2];
sx q[2];
rz(0.22656013) q[2];
rz(-1.4551) q[3];
sx q[3];
rz(-0.89613599) q[3];
sx q[3];
rz(-2.4494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9914472) q[0];
sx q[0];
rz(-1.8658072) q[0];
sx q[0];
rz(0.62514296) q[0];
rz(-1.9236247) q[1];
sx q[1];
rz(-0.70911276) q[1];
sx q[1];
rz(1.9295173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8344515) q[0];
sx q[0];
rz(-1.092415) q[0];
sx q[0];
rz(0.68104736) q[0];
x q[1];
rz(-1.3122968) q[2];
sx q[2];
rz(-0.91583672) q[2];
sx q[2];
rz(2.0531246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.068262488) q[1];
sx q[1];
rz(-2.0113328) q[1];
sx q[1];
rz(-0.063830062) q[1];
rz(-pi) q[2];
rz(1.8623705) q[3];
sx q[3];
rz(-1.8599038) q[3];
sx q[3];
rz(-1.5559199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43442279) q[2];
sx q[2];
rz(-1.3753128) q[2];
sx q[2];
rz(-2.0909069) q[2];
rz(0.55189842) q[3];
sx q[3];
rz(-0.69010186) q[3];
sx q[3];
rz(1.2845854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5156749) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(-2.3552409) q[1];
sx q[1];
rz(-2.4813589) q[1];
sx q[1];
rz(-2.9207041) q[1];
rz(-0.037842265) q[2];
sx q[2];
rz(-1.2076245) q[2];
sx q[2];
rz(-2.5943499) q[2];
rz(-3.0339387) q[3];
sx q[3];
rz(-0.68895491) q[3];
sx q[3];
rz(0.72910492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
