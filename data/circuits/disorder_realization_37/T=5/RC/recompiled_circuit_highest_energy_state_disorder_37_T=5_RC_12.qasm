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
rz(2.3647519) q[0];
sx q[0];
rz(-2.2651894) q[0];
sx q[0];
rz(-2.8887698) q[0];
rz(-1.9934935) q[1];
sx q[1];
rz(-0.91528457) q[1];
sx q[1];
rz(0.80438703) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0774121) q[0];
sx q[0];
rz(-1.8837116) q[0];
sx q[0];
rz(1.1699647) q[0];
rz(2.6635799) q[2];
sx q[2];
rz(-0.4768663) q[2];
sx q[2];
rz(3.06682) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6614395) q[1];
sx q[1];
rz(-2.6240908) q[1];
sx q[1];
rz(0.75019849) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7976296) q[3];
sx q[3];
rz(-1.4620145) q[3];
sx q[3];
rz(3.1358842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.048451) q[2];
sx q[2];
rz(-0.96258771) q[2];
sx q[2];
rz(-0.87240458) q[2];
rz(0.33556542) q[3];
sx q[3];
rz(-1.395697) q[3];
sx q[3];
rz(-0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63808477) q[0];
sx q[0];
rz(-2.8758949) q[0];
sx q[0];
rz(2.9124394) q[0];
rz(-0.19042641) q[1];
sx q[1];
rz(-1.3399905) q[1];
sx q[1];
rz(-0.71358877) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23799831) q[0];
sx q[0];
rz(-0.9146713) q[0];
sx q[0];
rz(-2.800992) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1515124) q[2];
sx q[2];
rz(-2.3186404) q[2];
sx q[2];
rz(-1.8320476) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.59871626) q[1];
sx q[1];
rz(-0.89756723) q[1];
sx q[1];
rz(-1.190541) q[1];
rz(0.17005597) q[3];
sx q[3];
rz(-2.3695558) q[3];
sx q[3];
rz(2.7626729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4552292) q[2];
sx q[2];
rz(-2.8296622) q[2];
sx q[2];
rz(0.4757821) q[2];
rz(-2.0717715) q[3];
sx q[3];
rz(-2.1333623) q[3];
sx q[3];
rz(0.76655918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0692724) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(1.2323761) q[0];
rz(-2.6909289) q[1];
sx q[1];
rz(-1.9602937) q[1];
sx q[1];
rz(2.252069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7256732) q[0];
sx q[0];
rz(-1.9751722) q[0];
sx q[0];
rz(-0.63979353) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3279104) q[2];
sx q[2];
rz(-0.70975862) q[2];
sx q[2];
rz(-0.26772945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3865859) q[1];
sx q[1];
rz(-2.4535547) q[1];
sx q[1];
rz(-1.2351456) q[1];
x q[2];
rz(-1.9185478) q[3];
sx q[3];
rz(-1.9257716) q[3];
sx q[3];
rz(-0.53670151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.687261) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(-1.2588151) q[2];
rz(1.9076294) q[3];
sx q[3];
rz(-2.0089269) q[3];
sx q[3];
rz(3.1335355) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4312129) q[0];
sx q[0];
rz(-0.90573913) q[0];
sx q[0];
rz(-1.4696962) q[0];
rz(-1.1471033) q[1];
sx q[1];
rz(-2.1397782) q[1];
sx q[1];
rz(-0.9009487) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44997901) q[0];
sx q[0];
rz(-1.8834582) q[0];
sx q[0];
rz(-1.6603966) q[0];
rz(-pi) q[1];
x q[1];
rz(2.35975) q[2];
sx q[2];
rz(-2.3580209) q[2];
sx q[2];
rz(1.9279566) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.683953) q[1];
sx q[1];
rz(-1.6648653) q[1];
sx q[1];
rz(2.4165618) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4086558) q[3];
sx q[3];
rz(-1.4519435) q[3];
sx q[3];
rz(-2.262923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49742302) q[2];
sx q[2];
rz(-1.7449417) q[2];
sx q[2];
rz(-2.4141342) q[2];
rz(2.4265031) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(-0.49811825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6775976) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(2.4053307) q[0];
rz(-0.62034208) q[1];
sx q[1];
rz(-2.5963929) q[1];
sx q[1];
rz(-1.6166519) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6432366) q[0];
sx q[0];
rz(-1.5806206) q[0];
sx q[0];
rz(3.0863484) q[0];
x q[1];
rz(-2.1651046) q[2];
sx q[2];
rz(-0.39820489) q[2];
sx q[2];
rz(-1.9081685) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.013629524) q[1];
sx q[1];
rz(-1.9444325) q[1];
sx q[1];
rz(-0.84898265) q[1];
x q[2];
rz(-0.78973887) q[3];
sx q[3];
rz(-2.2135255) q[3];
sx q[3];
rz(-1.5830884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5357431) q[2];
sx q[2];
rz(-0.77631408) q[2];
sx q[2];
rz(-2.1898451) q[2];
rz(-0.4176628) q[3];
sx q[3];
rz(-0.2499191) q[3];
sx q[3];
rz(-1.9865659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34053892) q[0];
sx q[0];
rz(-1.2507573) q[0];
sx q[0];
rz(2.387555) q[0];
rz(-2.2484089) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(2.975614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245859) q[0];
sx q[0];
rz(-0.69147325) q[0];
sx q[0];
rz(-2.9777479) q[0];
x q[1];
rz(2.6352623) q[2];
sx q[2];
rz(-1.4131318) q[2];
sx q[2];
rz(2.8540996) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9050203) q[1];
sx q[1];
rz(-0.81099866) q[1];
sx q[1];
rz(2.503997) q[1];
rz(-pi) q[2];
rz(-2.0319489) q[3];
sx q[3];
rz(-0.30546092) q[3];
sx q[3];
rz(-0.74126727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37504998) q[2];
sx q[2];
rz(-2.0039717) q[2];
sx q[2];
rz(3.1362015) q[2];
rz(3.052875) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(2.3458792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2583112) q[0];
sx q[0];
rz(-1.9323876) q[0];
sx q[0];
rz(-0.2051556) q[0];
rz(2.2329277) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(-0.044513449) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37986792) q[0];
sx q[0];
rz(-1.4141448) q[0];
sx q[0];
rz(-3.0101315) q[0];
x q[1];
rz(-1.8327863) q[2];
sx q[2];
rz(-2.9963608) q[2];
sx q[2];
rz(-2.0042208) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5959357) q[1];
sx q[1];
rz(-0.6164729) q[1];
sx q[1];
rz(1.193154) q[1];
rz(3.1188585) q[3];
sx q[3];
rz(-1.6904545) q[3];
sx q[3];
rz(2.5496063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8160416) q[2];
sx q[2];
rz(-2.6191235) q[2];
sx q[2];
rz(0.52118707) q[2];
rz(-2.9442287) q[3];
sx q[3];
rz(-1.3857931) q[3];
sx q[3];
rz(-2.638773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407783) q[0];
sx q[0];
rz(-0.1939119) q[0];
sx q[0];
rz(2.8863696) q[0];
rz(2.9995645) q[1];
sx q[1];
rz(-1.4165001) q[1];
sx q[1];
rz(-2.7878888) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3981374) q[0];
sx q[0];
rz(-2.6530735) q[0];
sx q[0];
rz(0.81870458) q[0];
rz(-pi) q[1];
rz(2.2204705) q[2];
sx q[2];
rz(-0.73474681) q[2];
sx q[2];
rz(-0.91780797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6981909) q[1];
sx q[1];
rz(-2.576722) q[1];
sx q[1];
rz(-2.2400212) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59732893) q[3];
sx q[3];
rz(-2.2976365) q[3];
sx q[3];
rz(-2.8387808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3069309) q[2];
sx q[2];
rz(-0.31535172) q[2];
sx q[2];
rz(1.6174512) q[2];
rz(2.2751685) q[3];
sx q[3];
rz(-1.4640936) q[3];
sx q[3];
rz(-0.58969897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6397112) q[0];
sx q[0];
rz(-2.9852133) q[0];
sx q[0];
rz(-1.3524652) q[0];
rz(1.9068708) q[1];
sx q[1];
rz(-0.23209485) q[1];
sx q[1];
rz(2.629705) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3054661) q[0];
sx q[0];
rz(-0.16634596) q[0];
sx q[0];
rz(1.8977099) q[0];
rz(0.16361605) q[2];
sx q[2];
rz(-2.3694042) q[2];
sx q[2];
rz(2.8761187) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7929188) q[1];
sx q[1];
rz(-1.9935384) q[1];
sx q[1];
rz(1.2482367) q[1];
rz(-1.0591577) q[3];
sx q[3];
rz(-2.4928204) q[3];
sx q[3];
rz(1.2973451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1341683) q[2];
sx q[2];
rz(-1.8210501) q[2];
sx q[2];
rz(-0.39719886) q[2];
rz(-2.126179) q[3];
sx q[3];
rz(-1.0011324) q[3];
sx q[3];
rz(-2.5889682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39887244) q[0];
sx q[0];
rz(-1.2074559) q[0];
sx q[0];
rz(0.43750986) q[0];
rz(-0.73879009) q[1];
sx q[1];
rz(-2.3414325) q[1];
sx q[1];
rz(-1.4791666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2302983) q[0];
sx q[0];
rz(-2.0574967) q[0];
sx q[0];
rz(-1.6337956) q[0];
x q[1];
rz(-1.9301038) q[2];
sx q[2];
rz(-1.6793161) q[2];
sx q[2];
rz(1.6720275) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65867546) q[1];
sx q[1];
rz(-1.5086996) q[1];
sx q[1];
rz(-3.0351522) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3607632) q[3];
sx q[3];
rz(-2.1159647) q[3];
sx q[3];
rz(0.050267537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2207569) q[2];
sx q[2];
rz(-1.3584542) q[2];
sx q[2];
rz(2.9662568) q[2];
rz(1.0026503) q[3];
sx q[3];
rz(-2.9084539) q[3];
sx q[3];
rz(1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47068448) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(0.15057527) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(0.038240959) q[2];
sx q[2];
rz(-0.19842166) q[2];
sx q[2];
rz(-1.7501065) q[2];
rz(0.2507052) q[3];
sx q[3];
rz(-1.3450932) q[3];
sx q[3];
rz(0.20592207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
