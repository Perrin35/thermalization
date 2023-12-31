OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(-1.0597205) q[0];
sx q[0];
rz(0.73097316) q[0];
rz(-1.5001186) q[1];
sx q[1];
rz(-2.1067696) q[1];
sx q[1];
rz(-2.1980481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9991) q[0];
sx q[0];
rz(-1.5038135) q[0];
sx q[0];
rz(-1.5357369) q[0];
rz(0.65096345) q[2];
sx q[2];
rz(-1.5032094) q[2];
sx q[2];
rz(-2.0219321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.01602068) q[1];
sx q[1];
rz(-1.9106094) q[1];
sx q[1];
rz(-1.1556975) q[1];
rz(2.5991873) q[3];
sx q[3];
rz(-2.7895658) q[3];
sx q[3];
rz(-2.5792518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.25508183) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(1.250766) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(-0.9799408) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.132906) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(-2.5426478) q[0];
rz(1.8006181) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(2.1751931) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67993977) q[0];
sx q[0];
rz(-1.1302117) q[0];
sx q[0];
rz(-0.64042129) q[0];
rz(-pi) q[1];
rz(0.62218372) q[2];
sx q[2];
rz(-1.6112279) q[2];
sx q[2];
rz(-1.7379023) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2981373) q[1];
sx q[1];
rz(-1.4569439) q[1];
sx q[1];
rz(1.4911806) q[1];
rz(-2.6838449) q[3];
sx q[3];
rz(-2.3380937) q[3];
sx q[3];
rz(2.1645434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0559343) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-0.30109626) q[2];
rz(-1.1931233) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0369204) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(-1.954129) q[0];
rz(-1.2359515) q[1];
sx q[1];
rz(-1.0373479) q[1];
sx q[1];
rz(-1.3175861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353969) q[0];
sx q[0];
rz(-0.72512308) q[0];
sx q[0];
rz(-0.32307415) q[0];
rz(-pi) q[1];
rz(-2.927711) q[2];
sx q[2];
rz(-1.3777395) q[2];
sx q[2];
rz(-2.4436827) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8220362) q[1];
sx q[1];
rz(-1.6720547) q[1];
sx q[1];
rz(-0.42145573) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0242689) q[3];
sx q[3];
rz(-2.085272) q[3];
sx q[3];
rz(-2.8770212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(-0.2066361) q[2];
rz(-0.7080428) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(-0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3669423) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(-2.2553717) q[0];
rz(-2.1318502) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(-1.2264235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3564295) q[0];
sx q[0];
rz(-0.7298846) q[0];
sx q[0];
rz(1.4714144) q[0];
x q[1];
rz(-2.2573651) q[2];
sx q[2];
rz(-1.9366169) q[2];
sx q[2];
rz(-2.5522752) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46982161) q[1];
sx q[1];
rz(-1.7676395) q[1];
sx q[1];
rz(-0.21398869) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63663441) q[3];
sx q[3];
rz(-2.5599179) q[3];
sx q[3];
rz(-0.10891529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4975138) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(-0.99299661) q[2];
rz(-1.3126866) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(-1.1859878) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.490373) q[1];
sx q[1];
rz(2.5591154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5929778) q[0];
sx q[0];
rz(-0.60615221) q[0];
sx q[0];
rz(-1.6968326) q[0];
rz(-pi) q[1];
rz(1.3864473) q[2];
sx q[2];
rz(-0.68612387) q[2];
sx q[2];
rz(-2.9550936) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1479668) q[1];
sx q[1];
rz(-2.2851351) q[1];
sx q[1];
rz(1.5453969) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4062823) q[3];
sx q[3];
rz(-1.5092106) q[3];
sx q[3];
rz(2.6046034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65486583) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(-0.081710903) q[2];
rz(2.667526) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896647) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(3.1337877) q[0];
rz(1.4004978) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(-1.1046462) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7846851) q[0];
sx q[0];
rz(-2.4570358) q[0];
sx q[0];
rz(0.25951578) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8203037) q[2];
sx q[2];
rz(-0.66547223) q[2];
sx q[2];
rz(-1.9208391) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.066597477) q[1];
sx q[1];
rz(-1.0269594) q[1];
sx q[1];
rz(2.4165513) q[1];
rz(-pi) q[2];
x q[2];
rz(2.036318) q[3];
sx q[3];
rz(-1.421531) q[3];
sx q[3];
rz(-2.535459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2509987) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(-2.5863623) q[2];
rz(0.17272078) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.5531439) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(3.0798262) q[0];
rz(-2.899509) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(-1.1118836) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2568946) q[0];
sx q[0];
rz(-1.9848794) q[0];
sx q[0];
rz(-0.82722442) q[0];
x q[1];
rz(-1.493181) q[2];
sx q[2];
rz(-1.5845808) q[2];
sx q[2];
rz(0.26253653) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4278533) q[1];
sx q[1];
rz(-0.76875988) q[1];
sx q[1];
rz(-0.40366918) q[1];
x q[2];
rz(0.63515969) q[3];
sx q[3];
rz(-2.7413462) q[3];
sx q[3];
rz(-1.1349585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8093439) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(2.3279482) q[2];
rz(1.404473) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(-2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(-1.7161436) q[0];
rz(-1.6199934) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(0.63751784) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2256048) q[0];
sx q[0];
rz(-0.84945744) q[0];
sx q[0];
rz(-2.3774873) q[0];
x q[1];
rz(-1.2136739) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(1.1536319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1453104) q[1];
sx q[1];
rz(-1.9687555) q[1];
sx q[1];
rz(2.5064962) q[1];
rz(-1.7082801) q[3];
sx q[3];
rz(-1.004389) q[3];
sx q[3];
rz(-0.24368225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0104388) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(-1.1361702) q[2];
rz(1.4853959) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(1.2095399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6256325) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(0.67614722) q[0];
rz(0.82538429) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(2.1264145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43468201) q[0];
sx q[0];
rz(-2.7043531) q[0];
sx q[0];
rz(0.31243639) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19998156) q[2];
sx q[2];
rz(-2.4790384) q[2];
sx q[2];
rz(2.8560864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0820513) q[1];
sx q[1];
rz(-2.6350628) q[1];
sx q[1];
rz(-0.7854714) q[1];
x q[2];
rz(2.7916662) q[3];
sx q[3];
rz(-0.79777989) q[3];
sx q[3];
rz(-1.5087023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72145808) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(0.96735111) q[2];
rz(-1.5445276) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(-2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6185146) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(3.0850947) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.2697376) q[1];
sx q[1];
rz(1.7369695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0987941) q[0];
sx q[0];
rz(-2.7569175) q[0];
sx q[0];
rz(-1.9210451) q[0];
rz(-pi) q[1];
rz(-2.769906) q[2];
sx q[2];
rz(-2.8166397) q[2];
sx q[2];
rz(1.5124958) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2056634) q[1];
sx q[1];
rz(-2.231039) q[1];
sx q[1];
rz(0.6026938) q[1];
x q[2];
rz(-0.057617188) q[3];
sx q[3];
rz(-2.9908097) q[3];
sx q[3];
rz(-2.5654716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29356062) q[2];
sx q[2];
rz(-2.3876987) q[2];
sx q[2];
rz(-2.8354697) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(2.117363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33481471) q[0];
sx q[0];
rz(-2.0756742) q[0];
sx q[0];
rz(-2.6609127) q[0];
rz(1.1595935) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(1.2724614) q[2];
sx q[2];
rz(-1.3186426) q[2];
sx q[2];
rz(-1.1124055) q[2];
rz(2.0648099) q[3];
sx q[3];
rz(-1.2772588) q[3];
sx q[3];
rz(0.96989934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
