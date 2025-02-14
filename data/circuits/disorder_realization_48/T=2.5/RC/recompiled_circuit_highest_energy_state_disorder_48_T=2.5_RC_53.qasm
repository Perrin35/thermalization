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
rz(2.7895522) q[0];
sx q[0];
rz(3.9665251) q[0];
sx q[0];
rz(12.031603) q[0];
rz(0.98400247) q[1];
sx q[1];
rz(0.50552955) q[1];
sx q[1];
rz(11.304392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.004744) q[0];
sx q[0];
rz(-1.8333577) q[0];
sx q[0];
rz(-2.5181818) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1970799) q[2];
sx q[2];
rz(-1.1188095) q[2];
sx q[2];
rz(-2.1747957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3024276) q[1];
sx q[1];
rz(-1.5189956) q[1];
sx q[1];
rz(-2.5940345) q[1];
x q[2];
rz(-2.7522911) q[3];
sx q[3];
rz(-1.3817781) q[3];
sx q[3];
rz(-1.0396599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24903211) q[2];
sx q[2];
rz(-0.51237115) q[2];
sx q[2];
rz(1.4830291) q[2];
rz(-2.856971) q[3];
sx q[3];
rz(-2.5183545) q[3];
sx q[3];
rz(2.0774138) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677218) q[0];
sx q[0];
rz(-2.4148648) q[0];
sx q[0];
rz(0.04100767) q[0];
rz(1.2076591) q[1];
sx q[1];
rz(-2.9137847) q[1];
sx q[1];
rz(-1.4606732) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8475094) q[0];
sx q[0];
rz(-0.14862157) q[0];
sx q[0];
rz(2.9411208) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3744205) q[2];
sx q[2];
rz(-1.1572654) q[2];
sx q[2];
rz(2.3656379) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41439287) q[1];
sx q[1];
rz(-1.314394) q[1];
sx q[1];
rz(3.1016769) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8051265) q[3];
sx q[3];
rz(-2.4357492) q[3];
sx q[3];
rz(0.88283759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66551208) q[2];
sx q[2];
rz(-1.4613232) q[2];
sx q[2];
rz(1.7512789) q[2];
rz(0.0811854) q[3];
sx q[3];
rz(-2.472671) q[3];
sx q[3];
rz(1.4359052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3944893) q[0];
sx q[0];
rz(-0.012270027) q[0];
sx q[0];
rz(2.5315206) q[0];
rz(-1.3129781) q[1];
sx q[1];
rz(-0.6956296) q[1];
sx q[1];
rz(-2.5909766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074215502) q[0];
sx q[0];
rz(-0.91892892) q[0];
sx q[0];
rz(-1.1135191) q[0];
rz(-pi) q[1];
rz(0.32260334) q[2];
sx q[2];
rz(-1.3122953) q[2];
sx q[2];
rz(-0.58835627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0455835) q[1];
sx q[1];
rz(-0.35638816) q[1];
sx q[1];
rz(1.0046887) q[1];
rz(-1.4161701) q[3];
sx q[3];
rz(-0.75115132) q[3];
sx q[3];
rz(0.42645833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6637471) q[2];
sx q[2];
rz(-0.17243324) q[2];
sx q[2];
rz(-2.0752068) q[2];
rz(1.9829228) q[3];
sx q[3];
rz(-1.7585157) q[3];
sx q[3];
rz(0.042044736) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31024194) q[0];
sx q[0];
rz(-0.61315918) q[0];
sx q[0];
rz(-0.9170652) q[0];
rz(2.4861368) q[1];
sx q[1];
rz(-0.90924811) q[1];
sx q[1];
rz(2.7520032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.675466) q[0];
sx q[0];
rz(-1.4891169) q[0];
sx q[0];
rz(1.7577175) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55928237) q[2];
sx q[2];
rz(-0.96508316) q[2];
sx q[2];
rz(-2.9756211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45976156) q[1];
sx q[1];
rz(-1.9719405) q[1];
sx q[1];
rz(0.031706867) q[1];
rz(-2.5960931) q[3];
sx q[3];
rz(-2.0974468) q[3];
sx q[3];
rz(-0.57034501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3098844) q[2];
sx q[2];
rz(-2.1990621) q[2];
sx q[2];
rz(1.7892828) q[2];
rz(2.3934707) q[3];
sx q[3];
rz(-1.8375405) q[3];
sx q[3];
rz(1.6149394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378976) q[0];
sx q[0];
rz(-2.5717773) q[0];
sx q[0];
rz(-1.8001528) q[0];
rz(2.1603284) q[1];
sx q[1];
rz(-1.2813247) q[1];
sx q[1];
rz(-2.431638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2316596) q[0];
sx q[0];
rz(-1.8965479) q[0];
sx q[0];
rz(1.5506106) q[0];
x q[1];
rz(2.9030062) q[2];
sx q[2];
rz(-1.0077395) q[2];
sx q[2];
rz(-2.4679139) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.73588307) q[1];
sx q[1];
rz(-0.27389474) q[1];
sx q[1];
rz(0.89400684) q[1];
rz(-pi) q[2];
rz(-0.59779915) q[3];
sx q[3];
rz(-1.4922659) q[3];
sx q[3];
rz(0.30818916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95825163) q[2];
sx q[2];
rz(-0.7603344) q[2];
sx q[2];
rz(0.90671268) q[2];
rz(0.18236154) q[3];
sx q[3];
rz(-1.4401108) q[3];
sx q[3];
rz(2.2831634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7957423) q[0];
sx q[0];
rz(-0.95473552) q[0];
sx q[0];
rz(2.9398651) q[0];
rz(-1.3564823) q[1];
sx q[1];
rz(-2.4701665) q[1];
sx q[1];
rz(1.9904402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1930453) q[0];
sx q[0];
rz(-2.6631115) q[0];
sx q[0];
rz(1.1043666) q[0];
rz(-pi) q[1];
rz(-1.1100805) q[2];
sx q[2];
rz(-2.1270129) q[2];
sx q[2];
rz(-1.0413025) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56111873) q[1];
sx q[1];
rz(-2.4137133) q[1];
sx q[1];
rz(2.5125458) q[1];
x q[2];
rz(-2.6309507) q[3];
sx q[3];
rz(-1.770442) q[3];
sx q[3];
rz(-0.24195237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1543033) q[2];
sx q[2];
rz(-0.58649784) q[2];
sx q[2];
rz(2.7941373) q[2];
rz(0.82184982) q[3];
sx q[3];
rz(-1.3653267) q[3];
sx q[3];
rz(-1.52389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3407985) q[0];
sx q[0];
rz(-0.64685416) q[0];
sx q[0];
rz(2.0183753) q[0];
rz(-0.42876354) q[1];
sx q[1];
rz(-2.2518497) q[1];
sx q[1];
rz(-0.14895983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70388795) q[0];
sx q[0];
rz(-1.0539248) q[0];
sx q[0];
rz(-0.025512841) q[0];
rz(-0.5908509) q[2];
sx q[2];
rz(-0.35195165) q[2];
sx q[2];
rz(0.9443493) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3330482) q[1];
sx q[1];
rz(-0.77189779) q[1];
sx q[1];
rz(0.21265642) q[1];
rz(-pi) q[2];
rz(1.4895053) q[3];
sx q[3];
rz(-0.40697655) q[3];
sx q[3];
rz(-1.5291884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1296156) q[2];
sx q[2];
rz(-0.35795438) q[2];
sx q[2];
rz(2.820106) q[2];
rz(1.1164411) q[3];
sx q[3];
rz(-2.2827086) q[3];
sx q[3];
rz(-1.727625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1534934) q[0];
sx q[0];
rz(-1.3445925) q[0];
sx q[0];
rz(3.1186812) q[0];
rz(1.5921536) q[1];
sx q[1];
rz(-1.0433082) q[1];
sx q[1];
rz(-1.2124088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1003636) q[0];
sx q[0];
rz(-2.1080906) q[0];
sx q[0];
rz(0.30502086) q[0];
x q[1];
rz(1.171807) q[2];
sx q[2];
rz(-2.3697457) q[2];
sx q[2];
rz(0.10153025) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.089253292) q[1];
sx q[1];
rz(-1.8732548) q[1];
sx q[1];
rz(-0.41771981) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0453835) q[3];
sx q[3];
rz(-0.85959496) q[3];
sx q[3];
rz(-2.8195153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93078485) q[2];
sx q[2];
rz(-0.19889861) q[2];
sx q[2];
rz(2.1442294) q[2];
rz(1.8831683) q[3];
sx q[3];
rz(-1.9505898) q[3];
sx q[3];
rz(-1.9112126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1579943) q[0];
sx q[0];
rz(-0.65584922) q[0];
sx q[0];
rz(2.6672145) q[0];
rz(-0.55167088) q[1];
sx q[1];
rz(-2.3730979) q[1];
sx q[1];
rz(2.4840568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5112572) q[0];
sx q[0];
rz(-1.5310107) q[0];
sx q[0];
rz(-1.5471094) q[0];
x q[1];
rz(-1.9507061) q[2];
sx q[2];
rz(-0.43336855) q[2];
sx q[2];
rz(0.72335183) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8847981) q[1];
sx q[1];
rz(-1.9391283) q[1];
sx q[1];
rz(-3.1241259) q[1];
rz(-3.0364584) q[3];
sx q[3];
rz(-1.7145559) q[3];
sx q[3];
rz(1.2105699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23003301) q[2];
sx q[2];
rz(-0.40731373) q[2];
sx q[2];
rz(1.8302906) q[2];
rz(-0.71410549) q[3];
sx q[3];
rz(-1.2823391) q[3];
sx q[3];
rz(0.56984058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1353726) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(-2.5332992) q[0];
rz(-2.3237806) q[1];
sx q[1];
rz(-1.1759956) q[1];
sx q[1];
rz(0.83293319) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21446709) q[0];
sx q[0];
rz(-1.9916849) q[0];
sx q[0];
rz(1.6273999) q[0];
x q[1];
rz(-1.3141749) q[2];
sx q[2];
rz(-1.875281) q[2];
sx q[2];
rz(-1.2198795) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6842839) q[1];
sx q[1];
rz(-2.5446518) q[1];
sx q[1];
rz(0.2376016) q[1];
rz(-pi) q[2];
rz(-0.59304955) q[3];
sx q[3];
rz(-1.5967973) q[3];
sx q[3];
rz(-3.0700695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1828764) q[2];
sx q[2];
rz(-0.35280886) q[2];
sx q[2];
rz(0.58615169) q[2];
rz(-2.2385249) q[3];
sx q[3];
rz(-2.51913) q[3];
sx q[3];
rz(-3.0557475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2387954) q[0];
sx q[0];
rz(-1.1302523) q[0];
sx q[0];
rz(-0.94373066) q[0];
rz(-1.0974274) q[1];
sx q[1];
rz(-2.8891017) q[1];
sx q[1];
rz(-0.15695922) q[1];
rz(3.0022754) q[2];
sx q[2];
rz(-1.2323772) q[2];
sx q[2];
rz(2.7215794) q[2];
rz(2.1395352) q[3];
sx q[3];
rz(-0.78227038) q[3];
sx q[3];
rz(-3.0113358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
