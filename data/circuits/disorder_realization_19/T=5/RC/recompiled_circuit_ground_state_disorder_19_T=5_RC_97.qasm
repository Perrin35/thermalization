OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.26025298) q[0];
sx q[0];
rz(-2.2007359) q[0];
sx q[0];
rz(0.22766222) q[0];
rz(-2.8582299) q[1];
sx q[1];
rz(-0.41937399) q[1];
sx q[1];
rz(-1.8546606) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0193699) q[0];
sx q[0];
rz(-1.625018) q[0];
sx q[0];
rz(2.0383459) q[0];
x q[1];
rz(-0.46704328) q[2];
sx q[2];
rz(-1.1640295) q[2];
sx q[2];
rz(-1.9860991) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70482774) q[1];
sx q[1];
rz(-2.4840762) q[1];
sx q[1];
rz(2.2566811) q[1];
rz(-2.5023978) q[3];
sx q[3];
rz(-1.3167644) q[3];
sx q[3];
rz(-2.3237236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8925573) q[2];
sx q[2];
rz(-2.1302569) q[2];
sx q[2];
rz(-2.2300143) q[2];
rz(2.3859731) q[3];
sx q[3];
rz(-0.32114649) q[3];
sx q[3];
rz(-2.6452981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3590473) q[0];
sx q[0];
rz(-1.8308715) q[0];
sx q[0];
rz(2.0027335) q[0];
rz(2.5149939) q[1];
sx q[1];
rz(-2.7732924) q[1];
sx q[1];
rz(-2.0638594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9288348) q[0];
sx q[0];
rz(-0.99883119) q[0];
sx q[0];
rz(1.1642745) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35524345) q[2];
sx q[2];
rz(-1.1525103) q[2];
sx q[2];
rz(2.929941) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30451286) q[1];
sx q[1];
rz(-2.8690352) q[1];
sx q[1];
rz(2.8902092) q[1];
x q[2];
rz(2.0502362) q[3];
sx q[3];
rz(-0.78244996) q[3];
sx q[3];
rz(-1.2408011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44800147) q[2];
sx q[2];
rz(-2.3841264) q[2];
sx q[2];
rz(2.9288911) q[2];
rz(1.5564144) q[3];
sx q[3];
rz(-0.74376619) q[3];
sx q[3];
rz(1.0630382) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0001275) q[0];
sx q[0];
rz(-0.26697049) q[0];
sx q[0];
rz(2.2947327) q[0];
rz(-0.2521387) q[1];
sx q[1];
rz(-0.82770258) q[1];
sx q[1];
rz(0.65021461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0668117) q[0];
sx q[0];
rz(-1.2322786) q[0];
sx q[0];
rz(-1.6174497) q[0];
x q[1];
rz(-2.3450801) q[2];
sx q[2];
rz(-1.0625417) q[2];
sx q[2];
rz(2.4426393) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1260155) q[1];
sx q[1];
rz(-1.460005) q[1];
sx q[1];
rz(-1.8689824) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3350983) q[3];
sx q[3];
rz(-1.8915081) q[3];
sx q[3];
rz(-0.93076462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3113159) q[2];
sx q[2];
rz(-0.69861424) q[2];
sx q[2];
rz(-0.96381956) q[2];
rz(1.7802995) q[3];
sx q[3];
rz(-2.6908974) q[3];
sx q[3];
rz(3.0986339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2394102) q[0];
sx q[0];
rz(-0.22265156) q[0];
sx q[0];
rz(-2.1233001) q[0];
rz(2.2564383) q[1];
sx q[1];
rz(-2.8926909) q[1];
sx q[1];
rz(-0.28908602) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8967696) q[0];
sx q[0];
rz(-0.45216783) q[0];
sx q[0];
rz(0.74613476) q[0];
x q[1];
rz(-1.9924329) q[2];
sx q[2];
rz(-2.1648063) q[2];
sx q[2];
rz(-0.27911148) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4319181) q[1];
sx q[1];
rz(-0.54952114) q[1];
sx q[1];
rz(1.9220244) q[1];
rz(-0.21896514) q[3];
sx q[3];
rz(-1.4462399) q[3];
sx q[3];
rz(-1.7809465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79124147) q[2];
sx q[2];
rz(-1.505932) q[2];
sx q[2];
rz(1.852847) q[2];
rz(0.57981235) q[3];
sx q[3];
rz(-0.47477397) q[3];
sx q[3];
rz(0.60540664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30699214) q[0];
sx q[0];
rz(-0.51486105) q[0];
sx q[0];
rz(2.8173764) q[0];
rz(0.038837198) q[1];
sx q[1];
rz(-0.7380929) q[1];
sx q[1];
rz(-2.0447581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34183472) q[0];
sx q[0];
rz(-1.8506764) q[0];
sx q[0];
rz(-1.3169692) q[0];
rz(-pi) q[1];
rz(3.039198) q[2];
sx q[2];
rz(-2.1490344) q[2];
sx q[2];
rz(1.6440934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8406163) q[1];
sx q[1];
rz(-2.8368717) q[1];
sx q[1];
rz(-1.4304377) q[1];
x q[2];
rz(1.7750793) q[3];
sx q[3];
rz(-2.7133803) q[3];
sx q[3];
rz(-0.34303676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.066112) q[2];
sx q[2];
rz(-0.26618633) q[2];
sx q[2];
rz(0.038662635) q[2];
rz(1.4085116) q[3];
sx q[3];
rz(-1.6950636) q[3];
sx q[3];
rz(-0.2814289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6283145) q[0];
sx q[0];
rz(-2.3230041) q[0];
sx q[0];
rz(2.1224838) q[0];
rz(0.45711532) q[1];
sx q[1];
rz(-2.8725084) q[1];
sx q[1];
rz(-1.4310744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4265617) q[0];
sx q[0];
rz(-2.9390897) q[0];
sx q[0];
rz(-3.0103548) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1598187) q[2];
sx q[2];
rz(-1.7830689) q[2];
sx q[2];
rz(0.50657099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0960321) q[1];
sx q[1];
rz(-1.211173) q[1];
sx q[1];
rz(-0.093861967) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39737153) q[3];
sx q[3];
rz(-2.2512967) q[3];
sx q[3];
rz(-1.3272663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9223601) q[2];
sx q[2];
rz(-1.7722426) q[2];
sx q[2];
rz(-0.57979453) q[2];
rz(-2.842105) q[3];
sx q[3];
rz(-2.6087285) q[3];
sx q[3];
rz(-1.9129491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2806468) q[0];
sx q[0];
rz(-1.4806643) q[0];
sx q[0];
rz(-0.10375599) q[0];
rz(-2.5170028) q[1];
sx q[1];
rz(-0.92085212) q[1];
sx q[1];
rz(-1.7594899) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.493295) q[0];
sx q[0];
rz(-1.7837423) q[0];
sx q[0];
rz(2.2981724) q[0];
x q[1];
rz(2.8417009) q[2];
sx q[2];
rz(-1.6768528) q[2];
sx q[2];
rz(2.7404355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15006062) q[1];
sx q[1];
rz(-2.074291) q[1];
sx q[1];
rz(0.30725422) q[1];
rz(-0.66094605) q[3];
sx q[3];
rz(-1.7978284) q[3];
sx q[3];
rz(1.1477136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78911191) q[2];
sx q[2];
rz(-1.4752957) q[2];
sx q[2];
rz(-2.2769807) q[2];
rz(-0.23884808) q[3];
sx q[3];
rz(-2.3819203) q[3];
sx q[3];
rz(2.4281832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9637941) q[0];
sx q[0];
rz(-0.56121427) q[0];
sx q[0];
rz(-2.906565) q[0];
rz(2.4759953) q[1];
sx q[1];
rz(-1.1382297) q[1];
sx q[1];
rz(-1.6682909) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4427933) q[0];
sx q[0];
rz(-2.5449552) q[0];
sx q[0];
rz(1.9477316) q[0];
rz(-pi) q[1];
rz(-1.7081882) q[2];
sx q[2];
rz(-1.1072888) q[2];
sx q[2];
rz(2.3799294) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43117796) q[1];
sx q[1];
rz(-1.3873552) q[1];
sx q[1];
rz(-1.6448433) q[1];
rz(-pi) q[2];
rz(0.45110945) q[3];
sx q[3];
rz(-1.530297) q[3];
sx q[3];
rz(-0.19865741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9792446) q[2];
sx q[2];
rz(-2.9408216) q[2];
sx q[2];
rz(0.46094224) q[2];
rz(-1.0848684) q[3];
sx q[3];
rz(-2.3960787) q[3];
sx q[3];
rz(2.357024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0775065) q[0];
sx q[0];
rz(-0.53090799) q[0];
sx q[0];
rz(-2.532646) q[0];
rz(-2.5755836) q[1];
sx q[1];
rz(-1.1606263) q[1];
sx q[1];
rz(-1.2014679) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58415675) q[0];
sx q[0];
rz(-3.0058001) q[0];
sx q[0];
rz(1.6312172) q[0];
x q[1];
rz(-1.0827093) q[2];
sx q[2];
rz(-1.7556453) q[2];
sx q[2];
rz(-2.3019997) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7215828) q[1];
sx q[1];
rz(-2.3370565) q[1];
sx q[1];
rz(-1.4415506) q[1];
rz(-pi) q[2];
rz(1.2668328) q[3];
sx q[3];
rz(-2.0018775) q[3];
sx q[3];
rz(1.6668741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0386049) q[2];
sx q[2];
rz(-1.0706341) q[2];
sx q[2];
rz(-1.9496244) q[2];
rz(-2.1206756) q[3];
sx q[3];
rz(-0.20191419) q[3];
sx q[3];
rz(1.5405704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9182619) q[0];
sx q[0];
rz(-0.20063618) q[0];
sx q[0];
rz(2.1740792) q[0];
rz(1.4406904) q[1];
sx q[1];
rz(-0.43165019) q[1];
sx q[1];
rz(-2.7593625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4995196) q[0];
sx q[0];
rz(-1.0040305) q[0];
sx q[0];
rz(0.21545179) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76127865) q[2];
sx q[2];
rz(-0.57028162) q[2];
sx q[2];
rz(-2.7340661) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8763242) q[1];
sx q[1];
rz(-2.3493715) q[1];
sx q[1];
rz(-1.4033615) q[1];
rz(-pi) q[2];
x q[2];
rz(0.02703826) q[3];
sx q[3];
rz(-2.595405) q[3];
sx q[3];
rz(2.9423713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.091577) q[2];
sx q[2];
rz(-0.86805582) q[2];
sx q[2];
rz(-0.78563219) q[2];
rz(0.026963726) q[3];
sx q[3];
rz(-1.5188768) q[3];
sx q[3];
rz(-0.54076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80326573) q[0];
sx q[0];
rz(-1.6852408) q[0];
sx q[0];
rz(2.2483873) q[0];
rz(-2.4772353) q[1];
sx q[1];
rz(-1.9021481) q[1];
sx q[1];
rz(-0.93217168) q[1];
rz(-1.2710089) q[2];
sx q[2];
rz(-2.2769417) q[2];
sx q[2];
rz(-1.4852038) q[2];
rz(-3.1278004) q[3];
sx q[3];
rz(-0.60473718) q[3];
sx q[3];
rz(-2.2040839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
