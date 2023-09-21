OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30480706) q[0];
sx q[0];
rz(-1.7903622) q[0];
sx q[0];
rz(0.027713393) q[0];
rz(2.2416441) q[2];
sx q[2];
rz(-0.64621011) q[2];
sx q[2];
rz(0.83713573) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5397415) q[1];
sx q[1];
rz(-1.1907693) q[1];
sx q[1];
rz(0.7976346) q[1];
rz(-pi) q[2];
rz(3.0580871) q[3];
sx q[3];
rz(-2.180035) q[3];
sx q[3];
rz(-3.1374251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41574079) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(-2.501781) q[2];
rz(-2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(-2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2663651) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(2.4282783) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(1.6289904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068802527) q[0];
sx q[0];
rz(-1.8913942) q[0];
sx q[0];
rz(0.94706236) q[0];
rz(-1.5741882) q[2];
sx q[2];
rz(-1.2870103) q[2];
sx q[2];
rz(1.1018745) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0257033) q[1];
sx q[1];
rz(-1.4374226) q[1];
sx q[1];
rz(-2.2504836) q[1];
x q[2];
rz(2.9019722) q[3];
sx q[3];
rz(-1.8488956) q[3];
sx q[3];
rz(0.77277377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(2.0837636) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(-3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2597044) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(-2.846068) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(0.74584109) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64622067) q[0];
sx q[0];
rz(-1.4797858) q[0];
sx q[0];
rz(1.792684) q[0];
x q[1];
rz(-0.27772851) q[2];
sx q[2];
rz(-0.54364294) q[2];
sx q[2];
rz(3.1090528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97400372) q[1];
sx q[1];
rz(-2.0473192) q[1];
sx q[1];
rz(-1.88489) q[1];
rz(0.48084534) q[3];
sx q[3];
rz(-2.6391609) q[3];
sx q[3];
rz(-1.06711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(2.4528465) q[2];
rz(0.050343242) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(-1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-2.922309) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(-3.0396089) q[0];
rz(-0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(2.8682958) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1647987) q[0];
sx q[0];
rz(-1.8935888) q[0];
sx q[0];
rz(-0.093683634) q[0];
rz(-2.1347413) q[2];
sx q[2];
rz(-2.1401005) q[2];
sx q[2];
rz(0.068892613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.30444333) q[1];
sx q[1];
rz(-1.7177561) q[1];
sx q[1];
rz(-1.8069581) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3917543) q[3];
sx q[3];
rz(-0.92902196) q[3];
sx q[3];
rz(0.87953506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(-2.6376574) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(0.082745634) q[0];
rz(2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-0.98714978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3521096) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(0.91870086) q[0];
x q[1];
rz(3.1290595) q[2];
sx q[2];
rz(-2.4977411) q[2];
sx q[2];
rz(-0.24521337) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5226749) q[1];
sx q[1];
rz(-2.4213311) q[1];
sx q[1];
rz(-0.48536761) q[1];
rz(-pi) q[2];
rz(3.0994814) q[3];
sx q[3];
rz(-1.1782421) q[3];
sx q[3];
rz(2.7819355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(-3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.485065) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(0.791839) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56023843) q[0];
sx q[0];
rz(-1.8525774) q[0];
sx q[0];
rz(0.015334107) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2811786) q[2];
sx q[2];
rz(-0.98266232) q[2];
sx q[2];
rz(-0.16298018) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2710323) q[1];
sx q[1];
rz(-2.6189657) q[1];
sx q[1];
rz(1.51103) q[1];
rz(-0.91026129) q[3];
sx q[3];
rz(-0.22781867) q[3];
sx q[3];
rz(-1.6662625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8906158) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(-0.77077579) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(-3.0859257) q[0];
rz(2.9260013) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(0.0035704426) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63821793) q[0];
sx q[0];
rz(-0.57045454) q[0];
sx q[0];
rz(0.72871448) q[0];
rz(-pi) q[1];
rz(1.9622757) q[2];
sx q[2];
rz(-2.2841095) q[2];
sx q[2];
rz(-1.6187402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7230941) q[1];
sx q[1];
rz(-1.2730036) q[1];
sx q[1];
rz(0.94350068) q[1];
x q[2];
rz(1.3725029) q[3];
sx q[3];
rz(-1.6453504) q[3];
sx q[3];
rz(-1.9004746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59721649) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6253117) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(2.4628941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5355276) q[0];
sx q[0];
rz(-1.8341656) q[0];
sx q[0];
rz(-0.33959629) q[0];
x q[1];
rz(2.3085824) q[2];
sx q[2];
rz(-1.2150803) q[2];
sx q[2];
rz(2.9290111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24819599) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(0.072695331) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42550605) q[3];
sx q[3];
rz(-1.7153499) q[3];
sx q[3];
rz(-1.0021462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17710182) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(0.93723047) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98638242) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(2.4329176) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-0.52694595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5567112) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(0.35993872) q[0];
rz(-pi) q[1];
rz(-1.2277463) q[2];
sx q[2];
rz(-0.68708778) q[2];
sx q[2];
rz(0.81817852) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2745797) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(-0.32151476) q[1];
rz(2.5832289) q[3];
sx q[3];
rz(-2.6287968) q[3];
sx q[3];
rz(-2.5496897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-0.33716831) q[3];
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
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(-2.1561484) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-0.3607761) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8737313) q[0];
sx q[0];
rz(-2.972313) q[0];
sx q[0];
rz(-1.7196359) q[0];
rz(-pi) q[1];
rz(-0.12014328) q[2];
sx q[2];
rz(-1.3573682) q[2];
sx q[2];
rz(-3.0394768) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0204748) q[1];
sx q[1];
rz(-1.4282244) q[1];
sx q[1];
rz(-2.3974182) q[1];
x q[2];
rz(0.1694069) q[3];
sx q[3];
rz(-1.425404) q[3];
sx q[3];
rz(-1.3225079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(-1.7760989) q[2];
sx q[2];
rz(-1.6663972) q[2];
sx q[2];
rz(-1.321928) q[2];
rz(-2.4102224) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
