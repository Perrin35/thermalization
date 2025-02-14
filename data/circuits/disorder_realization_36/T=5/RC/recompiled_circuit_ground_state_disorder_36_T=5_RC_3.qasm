OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1898243) q[0];
sx q[0];
rz(-0.70280743) q[0];
sx q[0];
rz(-1.9195317) q[0];
rz(3.089978) q[1];
sx q[1];
rz(-1.6366704) q[1];
sx q[1];
rz(1.6610891) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9771271) q[0];
sx q[0];
rz(-1.4071305) q[0];
sx q[0];
rz(2.2360691) q[0];
x q[1];
rz(2.1406581) q[2];
sx q[2];
rz(-2.9409932) q[2];
sx q[2];
rz(0.72884411) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.192321) q[1];
sx q[1];
rz(-1.1484129) q[1];
sx q[1];
rz(-1.7622838) q[1];
rz(-pi) q[2];
rz(2.7056115) q[3];
sx q[3];
rz(-2.0153305) q[3];
sx q[3];
rz(-1.3541927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9606955) q[2];
sx q[2];
rz(-0.81520671) q[2];
sx q[2];
rz(-1.0308456) q[2];
rz(2.893462) q[3];
sx q[3];
rz(-1.7957567) q[3];
sx q[3];
rz(-2.8743751) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2613075) q[0];
sx q[0];
rz(-1.945865) q[0];
sx q[0];
rz(-1.1241359) q[0];
rz(1.0719489) q[1];
sx q[1];
rz(-1.5644466) q[1];
sx q[1];
rz(-0.42676485) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65703553) q[0];
sx q[0];
rz(-0.74321514) q[0];
sx q[0];
rz(1.5667746) q[0];
rz(-0.34554225) q[2];
sx q[2];
rz(-1.8176518) q[2];
sx q[2];
rz(-0.051356476) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3755835) q[1];
sx q[1];
rz(-1.462862) q[1];
sx q[1];
rz(2.4426779) q[1];
rz(1.1077085) q[3];
sx q[3];
rz(-0.42965382) q[3];
sx q[3];
rz(0.11518726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4830604) q[2];
sx q[2];
rz(-0.51653647) q[2];
sx q[2];
rz(-3.0413682) q[2];
rz(-1.077486) q[3];
sx q[3];
rz(-1.5651549) q[3];
sx q[3];
rz(-1.235435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60233068) q[0];
sx q[0];
rz(-2.1644008) q[0];
sx q[0];
rz(-1.4343028) q[0];
rz(2.3151248) q[1];
sx q[1];
rz(-1.4974599) q[1];
sx q[1];
rz(-0.43063393) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9481107) q[0];
sx q[0];
rz(-1.3010345) q[0];
sx q[0];
rz(-1.7408235) q[0];
rz(-pi) q[1];
rz(-0.13027066) q[2];
sx q[2];
rz(-1.300087) q[2];
sx q[2];
rz(-0.66616466) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0525749) q[1];
sx q[1];
rz(-0.7886501) q[1];
sx q[1];
rz(-0.99747212) q[1];
rz(-0.96453676) q[3];
sx q[3];
rz(-2.5068847) q[3];
sx q[3];
rz(0.026594435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87381252) q[2];
sx q[2];
rz(-1.1649106) q[2];
sx q[2];
rz(2.5679892) q[2];
rz(-2.7576647) q[3];
sx q[3];
rz(-1.826518) q[3];
sx q[3];
rz(0.65281502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4189932) q[0];
sx q[0];
rz(-0.79177952) q[0];
sx q[0];
rz(2.0844039) q[0];
rz(0.81622299) q[1];
sx q[1];
rz(-0.99333251) q[1];
sx q[1];
rz(-2.9073471) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62158305) q[0];
sx q[0];
rz(-0.43061965) q[0];
sx q[0];
rz(-1.8355811) q[0];
rz(-pi) q[1];
rz(2.7170638) q[2];
sx q[2];
rz(-2.0574951) q[2];
sx q[2];
rz(-2.8667712) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.034742753) q[1];
sx q[1];
rz(-2.2487246) q[1];
sx q[1];
rz(2.2449916) q[1];
rz(1.0660421) q[3];
sx q[3];
rz(-1.1467993) q[3];
sx q[3];
rz(0.056886176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8546042) q[2];
sx q[2];
rz(-1.8082666) q[2];
sx q[2];
rz(-0.90281478) q[2];
rz(1.3826987) q[3];
sx q[3];
rz(-2.8185676) q[3];
sx q[3];
rz(-0.84669101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46120241) q[0];
sx q[0];
rz(-1.8620551) q[0];
sx q[0];
rz(0.70988208) q[0];
rz(-0.4711802) q[1];
sx q[1];
rz(-1.2268365) q[1];
sx q[1];
rz(-0.064362854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36633766) q[0];
sx q[0];
rz(-1.4357872) q[0];
sx q[0];
rz(-1.0827716) q[0];
x q[1];
rz(-2.9079403) q[2];
sx q[2];
rz(-1.086768) q[2];
sx q[2];
rz(2.7357227) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0546206) q[1];
sx q[1];
rz(-0.68212648) q[1];
sx q[1];
rz(-2.4404018) q[1];
rz(-pi) q[2];
rz(-1.5803171) q[3];
sx q[3];
rz(-1.4839982) q[3];
sx q[3];
rz(2.1824126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71089661) q[2];
sx q[2];
rz(-0.5527834) q[2];
sx q[2];
rz(-0.95332471) q[2];
rz(-2.583875) q[3];
sx q[3];
rz(-1.6744813) q[3];
sx q[3];
rz(-0.52887708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3144658) q[0];
sx q[0];
rz(-2.0966661) q[0];
sx q[0];
rz(-0.99660981) q[0];
rz(3.0820471) q[1];
sx q[1];
rz(-0.98809067) q[1];
sx q[1];
rz(-1.3522118) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72783453) q[0];
sx q[0];
rz(-1.0737077) q[0];
sx q[0];
rz(-0.97151206) q[0];
rz(2.8037595) q[2];
sx q[2];
rz(-0.80070904) q[2];
sx q[2];
rz(1.5986625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88207605) q[1];
sx q[1];
rz(-1.2959555) q[1];
sx q[1];
rz(-2.8311353) q[1];
x q[2];
rz(2.2289946) q[3];
sx q[3];
rz(-1.8479947) q[3];
sx q[3];
rz(2.4656221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.754564) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(0.55956364) q[2];
rz(-0.83545056) q[3];
sx q[3];
rz(-2.7159034) q[3];
sx q[3];
rz(1.3028418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.68592042) q[0];
sx q[0];
rz(-2.4113825) q[0];
sx q[0];
rz(-2.7749104) q[0];
rz(2.2600251) q[1];
sx q[1];
rz(-0.63789788) q[1];
sx q[1];
rz(-1.7104023) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2685705) q[0];
sx q[0];
rz(-2.2013469) q[0];
sx q[0];
rz(-1.0824758) q[0];
x q[1];
rz(-0.41421271) q[2];
sx q[2];
rz(-2.4117232) q[2];
sx q[2];
rz(3.102906) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.52932731) q[1];
sx q[1];
rz(-3.021486) q[1];
sx q[1];
rz(2.8976151) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4834385) q[3];
sx q[3];
rz(-1.7580716) q[3];
sx q[3];
rz(-1.8204456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7876579) q[2];
sx q[2];
rz(-0.5160318) q[2];
sx q[2];
rz(-2.3616974) q[2];
rz(2.6025313) q[3];
sx q[3];
rz(-1.5122248) q[3];
sx q[3];
rz(0.61162925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014538177) q[0];
sx q[0];
rz(-2.4712565) q[0];
sx q[0];
rz(0.77914733) q[0];
rz(0.51891333) q[1];
sx q[1];
rz(-1.2823558) q[1];
sx q[1];
rz(2.7386477) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084672734) q[0];
sx q[0];
rz(-0.92065996) q[0];
sx q[0];
rz(0.60303243) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8267118) q[2];
sx q[2];
rz(-2.0337031) q[2];
sx q[2];
rz(2.7969282) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0190474) q[1];
sx q[1];
rz(-1.9662939) q[1];
sx q[1];
rz(2.4707153) q[1];
x q[2];
rz(-0.99767606) q[3];
sx q[3];
rz(-1.3460288) q[3];
sx q[3];
rz(-2.5483957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0275823) q[2];
sx q[2];
rz(-1.1677914) q[2];
sx q[2];
rz(2.0243417) q[2];
rz(2.4452325) q[3];
sx q[3];
rz(-1.1098692) q[3];
sx q[3];
rz(1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47525147) q[0];
sx q[0];
rz(-0.23463686) q[0];
sx q[0];
rz(2.7269205) q[0];
rz(-2.0798202) q[1];
sx q[1];
rz(-2.2551408) q[1];
sx q[1];
rz(-0.83962238) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0939621) q[0];
sx q[0];
rz(-2.5174401) q[0];
sx q[0];
rz(2.312129) q[0];
x q[1];
rz(-0.7202486) q[2];
sx q[2];
rz(-2.0145085) q[2];
sx q[2];
rz(0.72138471) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66218195) q[1];
sx q[1];
rz(-2.4644797) q[1];
sx q[1];
rz(2.9138395) q[1];
rz(-pi) q[2];
rz(2.785061) q[3];
sx q[3];
rz(-1.339669) q[3];
sx q[3];
rz(0.69510776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35905579) q[2];
sx q[2];
rz(-1.7724719) q[2];
sx q[2];
rz(-3.0699406) q[2];
rz(-0.045529384) q[3];
sx q[3];
rz(-1.1979016) q[3];
sx q[3];
rz(-2.6825452) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5422106) q[0];
sx q[0];
rz(-2.5686503) q[0];
sx q[0];
rz(1.0546767) q[0];
rz(-0.032546267) q[1];
sx q[1];
rz(-2.1701505) q[1];
sx q[1];
rz(0.5221101) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3536026) q[0];
sx q[0];
rz(-0.66594571) q[0];
sx q[0];
rz(0.5086201) q[0];
x q[1];
rz(-0.73660518) q[2];
sx q[2];
rz(-2.0275368) q[2];
sx q[2];
rz(2.0219986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0070399) q[1];
sx q[1];
rz(-1.947548) q[1];
sx q[1];
rz(-1.9794078) q[1];
rz(0.22618146) q[3];
sx q[3];
rz(-1.0345248) q[3];
sx q[3];
rz(0.29335833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1471499) q[2];
sx q[2];
rz(-0.20211896) q[2];
sx q[2];
rz(-1.7333376) q[2];
rz(-1.555892) q[3];
sx q[3];
rz(-0.2318016) q[3];
sx q[3];
rz(-0.91109341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0178575) q[0];
sx q[0];
rz(-1.8386848) q[0];
sx q[0];
rz(0.89717502) q[0];
rz(-2.1625715) q[1];
sx q[1];
rz(-1.1430102) q[1];
sx q[1];
rz(2.7390726) q[1];
rz(2.0918905) q[2];
sx q[2];
rz(-2.7626531) q[2];
sx q[2];
rz(0.74041453) q[2];
rz(-2.4784563) q[3];
sx q[3];
rz(-1.7189773) q[3];
sx q[3];
rz(1.1130352) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
