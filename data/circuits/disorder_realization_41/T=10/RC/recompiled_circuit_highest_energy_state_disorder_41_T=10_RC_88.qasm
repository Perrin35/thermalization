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
rz(-2.8180985) q[0];
sx q[0];
rz(-2.9383724) q[0];
sx q[0];
rz(2.8414371) q[0];
rz(-0.35429859) q[1];
sx q[1];
rz(4.3759182) q[1];
sx q[1];
rz(10.195923) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010492652) q[0];
sx q[0];
rz(-0.26828921) q[0];
sx q[0];
rz(-2.3649251) q[0];
x q[1];
rz(-2.6080898) q[2];
sx q[2];
rz(-1.9868104) q[2];
sx q[2];
rz(-2.5737284) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59471666) q[1];
sx q[1];
rz(-1.1961814) q[1];
sx q[1];
rz(-2.0105816) q[1];
rz(-pi) q[2];
rz(-1.7158137) q[3];
sx q[3];
rz(-0.43802625) q[3];
sx q[3];
rz(-1.6159542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1511718) q[2];
sx q[2];
rz(-2.6800818) q[2];
sx q[2];
rz(1.0837466) q[2];
rz(-0.57461965) q[3];
sx q[3];
rz(-1.5950404) q[3];
sx q[3];
rz(0.77742022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3144749) q[0];
sx q[0];
rz(-0.62905335) q[0];
sx q[0];
rz(2.7431059) q[0];
rz(1.1442979) q[1];
sx q[1];
rz(-2.5356348) q[1];
sx q[1];
rz(0.95432895) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9480316) q[0];
sx q[0];
rz(-1.5662417) q[0];
sx q[0];
rz(-1.5713552) q[0];
rz(-3.047557) q[2];
sx q[2];
rz(-0.91528972) q[2];
sx q[2];
rz(1.2588225) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9056503) q[1];
sx q[1];
rz(-0.62726089) q[1];
sx q[1];
rz(-1.1494517) q[1];
x q[2];
rz(-2.6000836) q[3];
sx q[3];
rz(-1.1605822) q[3];
sx q[3];
rz(2.4334986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2362471) q[2];
sx q[2];
rz(-0.50062847) q[2];
sx q[2];
rz(3.0312669) q[2];
rz(-2.3162383) q[3];
sx q[3];
rz(-0.76485157) q[3];
sx q[3];
rz(-0.74025214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1684619) q[0];
sx q[0];
rz(-0.2178807) q[0];
sx q[0];
rz(-0.88200283) q[0];
rz(1.0285671) q[1];
sx q[1];
rz(-0.55362892) q[1];
sx q[1];
rz(2.2291768) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9498429) q[0];
sx q[0];
rz(-2.0217359) q[0];
sx q[0];
rz(-0.080206932) q[0];
rz(-pi) q[1];
rz(-1.4503612) q[2];
sx q[2];
rz(-0.85589534) q[2];
sx q[2];
rz(0.038106136) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0076651) q[1];
sx q[1];
rz(-1.0093736) q[1];
sx q[1];
rz(-0.39427653) q[1];
x q[2];
rz(1.1242718) q[3];
sx q[3];
rz(-1.339661) q[3];
sx q[3];
rz(2.6416306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.87938157) q[2];
sx q[2];
rz(-0.24803455) q[2];
sx q[2];
rz(-0.81825078) q[2];
rz(0.36000559) q[3];
sx q[3];
rz(-1.2986978) q[3];
sx q[3];
rz(-0.037809614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4358732) q[0];
sx q[0];
rz(-2.191045) q[0];
sx q[0];
rz(-1.9807504) q[0];
rz(0.97087651) q[1];
sx q[1];
rz(-0.91578805) q[1];
sx q[1];
rz(-0.11347778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.023983) q[0];
sx q[0];
rz(-1.1679839) q[0];
sx q[0];
rz(-0.65448032) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6062097) q[2];
sx q[2];
rz(-1.5839911) q[2];
sx q[2];
rz(-2.925417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7474688) q[1];
sx q[1];
rz(-1.040084) q[1];
sx q[1];
rz(2.3783724) q[1];
rz(-pi) q[2];
rz(0.63469074) q[3];
sx q[3];
rz(-1.029976) q[3];
sx q[3];
rz(1.5405003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11161441) q[2];
sx q[2];
rz(-1.8054211) q[2];
sx q[2];
rz(-2.2335936) q[2];
rz(2.8152483) q[3];
sx q[3];
rz(-0.55485266) q[3];
sx q[3];
rz(-2.7197796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.570785) q[0];
sx q[0];
rz(-0.76199216) q[0];
sx q[0];
rz(-2.1245891) q[0];
rz(2.6550338) q[1];
sx q[1];
rz(-1.0252955) q[1];
sx q[1];
rz(-3.0466381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61750283) q[0];
sx q[0];
rz(-1.4727778) q[0];
sx q[0];
rz(1.4630415) q[0];
rz(-pi) q[1];
rz(2.3170065) q[2];
sx q[2];
rz(-1.4207625) q[2];
sx q[2];
rz(-2.4730446) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1639171) q[1];
sx q[1];
rz(-1.9174121) q[1];
sx q[1];
rz(1.279128) q[1];
x q[2];
rz(-1.2566936) q[3];
sx q[3];
rz(-1.2746547) q[3];
sx q[3];
rz(1.2552798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.93593705) q[2];
sx q[2];
rz(-0.63465261) q[2];
sx q[2];
rz(0.46979365) q[2];
rz(2.4433686) q[3];
sx q[3];
rz(-1.9152812) q[3];
sx q[3];
rz(0.32018143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6093589) q[0];
sx q[0];
rz(-2.4025752) q[0];
sx q[0];
rz(-0.21482378) q[0];
rz(-2.9048982) q[1];
sx q[1];
rz(-2.5645945) q[1];
sx q[1];
rz(2.7223041) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0146862) q[0];
sx q[0];
rz(-1.4923411) q[0];
sx q[0];
rz(-0.04562745) q[0];
x q[1];
rz(-2.0620286) q[2];
sx q[2];
rz(-0.59796792) q[2];
sx q[2];
rz(0.21737013) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46696906) q[1];
sx q[1];
rz(-1.7745275) q[1];
sx q[1];
rz(1.2342888) q[1];
x q[2];
rz(-0.71281302) q[3];
sx q[3];
rz(-1.7216923) q[3];
sx q[3];
rz(3.0217882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99570167) q[2];
sx q[2];
rz(-0.13549165) q[2];
sx q[2];
rz(-3.0502012) q[2];
rz(2.7898096) q[3];
sx q[3];
rz(-2.3942949) q[3];
sx q[3];
rz(-2.6763797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9212937) q[0];
sx q[0];
rz(-1.969368) q[0];
sx q[0];
rz(-1.9660796) q[0];
rz(-0.57918864) q[1];
sx q[1];
rz(-0.80668956) q[1];
sx q[1];
rz(-2.9520292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2611147) q[0];
sx q[0];
rz(-1.3700598) q[0];
sx q[0];
rz(1.6693382) q[0];
x q[1];
rz(-1.4681508) q[2];
sx q[2];
rz(-0.70404875) q[2];
sx q[2];
rz(-1.8542022) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2793365) q[1];
sx q[1];
rz(-1.5784043) q[1];
sx q[1];
rz(-0.81806446) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0814384) q[3];
sx q[3];
rz(-1.0511479) q[3];
sx q[3];
rz(-2.1148256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2438948) q[2];
sx q[2];
rz(-1.8420668) q[2];
sx q[2];
rz(0.50590903) q[2];
rz(0.77887744) q[3];
sx q[3];
rz(-0.39611045) q[3];
sx q[3];
rz(1.8469384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8706354) q[0];
sx q[0];
rz(-2.0896235) q[0];
sx q[0];
rz(1.9860995) q[0];
rz(-0.0011860154) q[1];
sx q[1];
rz(-2.3540034) q[1];
sx q[1];
rz(0.40518618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14286539) q[0];
sx q[0];
rz(-2.0901276) q[0];
sx q[0];
rz(-2.1148391) q[0];
rz(1.1616906) q[2];
sx q[2];
rz(-1.7120581) q[2];
sx q[2];
rz(-0.90785949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0464814) q[1];
sx q[1];
rz(-0.80911923) q[1];
sx q[1];
rz(-3.0363096) q[1];
x q[2];
rz(-2.6514945) q[3];
sx q[3];
rz(-0.14486545) q[3];
sx q[3];
rz(-2.9165096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81795168) q[2];
sx q[2];
rz(-1.9743071) q[2];
sx q[2];
rz(-0.52158588) q[2];
rz(3.0449872) q[3];
sx q[3];
rz(-2.1229027) q[3];
sx q[3];
rz(0.49083138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7535962) q[0];
sx q[0];
rz(-0.8803519) q[0];
sx q[0];
rz(3.0269347) q[0];
rz(-1.8278587) q[1];
sx q[1];
rz(-1.6856245) q[1];
sx q[1];
rz(-3.0065261) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60914441) q[0];
sx q[0];
rz(-1.3125064) q[0];
sx q[0];
rz(-1.3468955) q[0];
rz(-1.8377322) q[2];
sx q[2];
rz(-1.1914754) q[2];
sx q[2];
rz(2.1960559) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5852063) q[1];
sx q[1];
rz(-1.7328246) q[1];
sx q[1];
rz(0.94604793) q[1];
rz(1.6879795) q[3];
sx q[3];
rz(-2.5796267) q[3];
sx q[3];
rz(-2.9148852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31585285) q[2];
sx q[2];
rz(-1.9261381) q[2];
sx q[2];
rz(-2.1960171) q[2];
rz(2.8816176) q[3];
sx q[3];
rz(-2.1187449) q[3];
sx q[3];
rz(2.0135349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.05397) q[0];
sx q[0];
rz(-2.051351) q[0];
sx q[0];
rz(-1.7386275) q[0];
rz(-2.2258017) q[1];
sx q[1];
rz(-2.5251838) q[1];
sx q[1];
rz(0.91555196) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2820194) q[0];
sx q[0];
rz(-1.440319) q[0];
sx q[0];
rz(3.0436712) q[0];
rz(-pi) q[1];
x q[1];
rz(2.948817) q[2];
sx q[2];
rz(-2.0844242) q[2];
sx q[2];
rz(2.7889268) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3183704) q[1];
sx q[1];
rz(-2.530066) q[1];
sx q[1];
rz(2.3182858) q[1];
rz(2.2804349) q[3];
sx q[3];
rz(-1.4415411) q[3];
sx q[3];
rz(-1.2354163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0182858) q[2];
sx q[2];
rz(-0.64299631) q[2];
sx q[2];
rz(2.7389738) q[2];
rz(-1.9649402) q[3];
sx q[3];
rz(-2.8549356) q[3];
sx q[3];
rz(-0.76796842) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9993512) q[0];
sx q[0];
rz(-0.9724697) q[0];
sx q[0];
rz(-1.019626) q[0];
rz(-2.453852) q[1];
sx q[1];
rz(-0.80324695) q[1];
sx q[1];
rz(-1.3245503) q[1];
rz(-2.4551212) q[2];
sx q[2];
rz(-2.1778637) q[2];
sx q[2];
rz(-1.7354497) q[2];
rz(1.8263893) q[3];
sx q[3];
rz(-2.1859313) q[3];
sx q[3];
rz(1.9060322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
