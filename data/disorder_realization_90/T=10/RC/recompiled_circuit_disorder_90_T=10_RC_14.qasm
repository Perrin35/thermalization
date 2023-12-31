OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(3.0545711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10843033) q[0];
sx q[0];
rz(-0.96956367) q[0];
sx q[0];
rz(0.45494672) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1331698) q[2];
sx q[2];
rz(-2.3699017) q[2];
sx q[2];
rz(-1.1222249) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.39057186) q[1];
sx q[1];
rz(-0.50288548) q[1];
sx q[1];
rz(-1.2807756) q[1];
rz(-pi) q[2];
rz(-1.617241) q[3];
sx q[3];
rz(-1.9825476) q[3];
sx q[3];
rz(2.9573033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0027851) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574629) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.194838) q[0];
rz(3.0589814) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(-0.00037489051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9276792) q[0];
sx q[0];
rz(-1.9880901) q[0];
sx q[0];
rz(-0.37950619) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92210251) q[2];
sx q[2];
rz(-0.32425913) q[2];
sx q[2];
rz(-2.1622554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2692711) q[1];
sx q[1];
rz(-1.806013) q[1];
sx q[1];
rz(-0.61551952) q[1];
rz(-3.0307816) q[3];
sx q[3];
rz(-2.6786945) q[3];
sx q[3];
rz(-3.0447931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80883819) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(-0.17671281) q[2];
rz(-0.79408944) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(-0.55364048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(-2.7368271) q[0];
rz(1.8602712) q[1];
sx q[1];
rz(-2.5679913) q[1];
sx q[1];
rz(-1.8331029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3151911) q[0];
sx q[0];
rz(-2.0668525) q[0];
sx q[0];
rz(1.0415003) q[0];
x q[1];
rz(-1.8070418) q[2];
sx q[2];
rz(-1.0457195) q[2];
sx q[2];
rz(0.46307785) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1723459) q[1];
sx q[1];
rz(-1.8241276) q[1];
sx q[1];
rz(-2.2971056) q[1];
rz(-0.43867302) q[3];
sx q[3];
rz(-1.4891426) q[3];
sx q[3];
rz(0.7338394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(3.1211839) q[2];
rz(1.071788) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(0.91598696) q[0];
rz(-0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(1.0571009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14295828) q[0];
sx q[0];
rz(-2.7035473) q[0];
sx q[0];
rz(-0.060158718) q[0];
x q[1];
rz(-0.57025036) q[2];
sx q[2];
rz(-1.0409365) q[2];
sx q[2];
rz(2.3392764) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4790198) q[1];
sx q[1];
rz(-2.9111007) q[1];
sx q[1];
rz(-1.6307994) q[1];
rz(-0.75745849) q[3];
sx q[3];
rz(-1.4043024) q[3];
sx q[3];
rz(1.7050626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2477734) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(-2.7988953) q[2];
rz(1.6977067) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-0.56348339) q[0];
sx q[0];
rz(3.0551531) q[0];
rz(-1.3899639) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(-2.6729029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5993354) q[0];
sx q[0];
rz(-0.89878966) q[0];
sx q[0];
rz(-0.50962781) q[0];
x q[1];
rz(1.3895949) q[2];
sx q[2];
rz(-1.2672408) q[2];
sx q[2];
rz(-0.62523491) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9710755) q[1];
sx q[1];
rz(-1.7726388) q[1];
sx q[1];
rz(-1.7032196) q[1];
rz(-pi) q[2];
rz(-0.080805578) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9995352) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(-0.86432499) q[2];
rz(-0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(-2.8488081) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7140759) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(2.9445904) q[0];
rz(1.3621832) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(0.33624712) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2339736) q[0];
sx q[0];
rz(-0.95171463) q[0];
sx q[0];
rz(1.0236077) q[0];
rz(-2.3669846) q[2];
sx q[2];
rz(-0.75658549) q[2];
sx q[2];
rz(0.87385439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.058192249) q[1];
sx q[1];
rz(-2.5501049) q[1];
sx q[1];
rz(-3.0565492) q[1];
rz(-2.9365262) q[3];
sx q[3];
rz(-0.70687095) q[3];
sx q[3];
rz(-2.6170078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.53529915) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(0.84632787) q[2];
rz(-1.9019295) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(-0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.7476615) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(-0.25892648) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(2.0475725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8261995) q[0];
sx q[0];
rz(-2.0059735) q[0];
sx q[0];
rz(2.7639593) q[0];
x q[1];
rz(1.5161683) q[2];
sx q[2];
rz(-0.22287456) q[2];
sx q[2];
rz(-0.54578997) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8158405) q[1];
sx q[1];
rz(-3.0430803) q[1];
sx q[1];
rz(1.2835531) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21906994) q[3];
sx q[3];
rz(-0.1004569) q[3];
sx q[3];
rz(-0.62517525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(0.32067498) q[2];
rz(0.43618068) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(-2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748579) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(2.136769) q[0];
rz(-2.5514305) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(-0.80387962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3636712) q[0];
sx q[0];
rz(-1.5313308) q[0];
sx q[0];
rz(2.5556593) q[0];
rz(-pi) q[1];
rz(-0.55401037) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(-0.49288921) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1357321) q[1];
sx q[1];
rz(-1.5605968) q[1];
sx q[1];
rz(-1.6594396) q[1];
x q[2];
rz(1.6611093) q[3];
sx q[3];
rz(-2.2323425) q[3];
sx q[3];
rz(-2.6722398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33264318) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-1.0104898) q[2];
rz(0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(-0.4075152) q[0];
rz(0.28911668) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-0.75072748) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6347046) q[0];
sx q[0];
rz(-1.2872211) q[0];
sx q[0];
rz(-0.30446913) q[0];
rz(-pi) q[1];
rz(-0.27751343) q[2];
sx q[2];
rz(-2.8515184) q[2];
sx q[2];
rz(1.2558503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1474485) q[1];
sx q[1];
rz(-1.7732883) q[1];
sx q[1];
rz(1.8840428) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0575033) q[3];
sx q[3];
rz(-0.30232271) q[3];
sx q[3];
rz(1.38894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0728545) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(1.1661952) q[2];
rz(1.3646305) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5584548) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(1.3903842) q[0];
rz(-0.33069262) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(1.4601382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18171039) q[0];
sx q[0];
rz(-1.0707756) q[0];
sx q[0];
rz(2.659003) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1183364) q[2];
sx q[2];
rz(-2.0746982) q[2];
sx q[2];
rz(1.1018673) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8148988) q[1];
sx q[1];
rz(-1.9683451) q[1];
sx q[1];
rz(1.7705998) q[1];
x q[2];
rz(-0.71698935) q[3];
sx q[3];
rz(-0.90962142) q[3];
sx q[3];
rz(-0.17718525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7918487) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(-1.4617408) q[2];
rz(2.0215624) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-0.45564836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6256975) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.3810146) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(1.21576) q[2];
sx q[2];
rz(-2.3264865) q[2];
sx q[2];
rz(-2.9441499) q[2];
rz(1.1827042) q[3];
sx q[3];
rz(-2.2065065) q[3];
sx q[3];
rz(1.828215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
