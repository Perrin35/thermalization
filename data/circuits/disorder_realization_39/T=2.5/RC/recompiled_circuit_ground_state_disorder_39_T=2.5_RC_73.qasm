OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.429739) q[0];
sx q[0];
rz(-1.9965594) q[0];
sx q[0];
rz(-1.0166919) q[0];
rz(-3.9859803) q[1];
sx q[1];
rz(3.9581668) q[1];
sx q[1];
rz(8.858455) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9670463) q[0];
sx q[0];
rz(-2.5264347) q[0];
sx q[0];
rz(-1.7394203) q[0];
rz(-pi) q[1];
rz(0.077458642) q[2];
sx q[2];
rz(-1.768154) q[2];
sx q[2];
rz(-3.0826621) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29222088) q[1];
sx q[1];
rz(-1.148889) q[1];
sx q[1];
rz(-1.6590236) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8352691) q[3];
sx q[3];
rz(-0.60738436) q[3];
sx q[3];
rz(-0.67061433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5113968) q[2];
sx q[2];
rz(-2.2679195) q[2];
sx q[2];
rz(-2.013618) q[2];
rz(1.6389716) q[3];
sx q[3];
rz(-2.2962544) q[3];
sx q[3];
rz(-0.42475167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86654919) q[0];
sx q[0];
rz(-0.4568704) q[0];
sx q[0];
rz(3.1006815) q[0];
rz(-2.5919137) q[1];
sx q[1];
rz(-0.37893852) q[1];
sx q[1];
rz(-0.080370195) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5588829) q[0];
sx q[0];
rz(-2.2930995) q[0];
sx q[0];
rz(1.5651907) q[0];
rz(-pi) q[1];
rz(2.2533312) q[2];
sx q[2];
rz(-0.76768926) q[2];
sx q[2];
rz(-1.7861799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3894905) q[1];
sx q[1];
rz(-2.6536149) q[1];
sx q[1];
rz(1.6882903) q[1];
rz(-pi) q[2];
rz(1.188813) q[3];
sx q[3];
rz(-2.502591) q[3];
sx q[3];
rz(-1.5843713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4139159) q[2];
sx q[2];
rz(-0.84299403) q[2];
sx q[2];
rz(-0.28398871) q[2];
rz(1.5208288) q[3];
sx q[3];
rz(-1.5512543) q[3];
sx q[3];
rz(-0.76333299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82336998) q[0];
sx q[0];
rz(-2.9871873) q[0];
sx q[0];
rz(-0.41686091) q[0];
rz(2.2082224) q[1];
sx q[1];
rz(-0.88637543) q[1];
sx q[1];
rz(-2.0534168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4517072) q[0];
sx q[0];
rz(-2.5390115) q[0];
sx q[0];
rz(-1.129928) q[0];
x q[1];
rz(-2.5809885) q[2];
sx q[2];
rz(-1.6050771) q[2];
sx q[2];
rz(2.3115668) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29050204) q[1];
sx q[1];
rz(-2.7298285) q[1];
sx q[1];
rz(-0.74300933) q[1];
x q[2];
rz(2.9253694) q[3];
sx q[3];
rz(-0.29400405) q[3];
sx q[3];
rz(2.2764429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0205959) q[2];
sx q[2];
rz(-0.72046295) q[2];
sx q[2];
rz(2.8159115) q[2];
rz(0.39564141) q[3];
sx q[3];
rz(-0.49042693) q[3];
sx q[3];
rz(2.4237848) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2242103) q[0];
sx q[0];
rz(-2.2069187) q[0];
sx q[0];
rz(-0.14064661) q[0];
rz(-2.7463101) q[1];
sx q[1];
rz(-1.544516) q[1];
sx q[1];
rz(-0.43221727) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062586322) q[0];
sx q[0];
rz(-2.5012928) q[0];
sx q[0];
rz(-2.2878134) q[0];
rz(-pi) q[1];
rz(1.4163966) q[2];
sx q[2];
rz(-2.6302166) q[2];
sx q[2];
rz(-1.7768154) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85104698) q[1];
sx q[1];
rz(-1.6509027) q[1];
sx q[1];
rz(-2.8677031) q[1];
x q[2];
rz(-2.3234576) q[3];
sx q[3];
rz(-2.3664118) q[3];
sx q[3];
rz(1.1213746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.846659) q[2];
sx q[2];
rz(-1.1120956) q[2];
sx q[2];
rz(2.6726932) q[2];
rz(2.9851959) q[3];
sx q[3];
rz(-2.1726435) q[3];
sx q[3];
rz(-1.8739353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073931996) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(-0.78039783) q[0];
rz(2.228915) q[1];
sx q[1];
rz(-2.3527805) q[1];
sx q[1];
rz(1.7524293) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71043832) q[0];
sx q[0];
rz(-1.2200933) q[0];
sx q[0];
rz(1.2855269) q[0];
rz(-pi) q[1];
rz(1.7002128) q[2];
sx q[2];
rz(-0.99301087) q[2];
sx q[2];
rz(-0.31020626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.011574419) q[1];
sx q[1];
rz(-2.5352074) q[1];
sx q[1];
rz(2.9467086) q[1];
rz(-pi) q[2];
rz(-2.6776047) q[3];
sx q[3];
rz(-0.9404805) q[3];
sx q[3];
rz(0.88018552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7954365) q[2];
sx q[2];
rz(-2.4861591) q[2];
sx q[2];
rz(-2.9620841) q[2];
rz(1.6620212) q[3];
sx q[3];
rz(-0.69412762) q[3];
sx q[3];
rz(0.0050541335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805304) q[0];
sx q[0];
rz(-1.625165) q[0];
sx q[0];
rz(1.6777212) q[0];
rz(0.013710984) q[1];
sx q[1];
rz(-1.6746215) q[1];
sx q[1];
rz(0.94508583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4714519) q[0];
sx q[0];
rz(-1.9197122) q[0];
sx q[0];
rz(-0.80100153) q[0];
rz(-pi) q[1];
rz(-2.5588972) q[2];
sx q[2];
rz(-2.7740363) q[2];
sx q[2];
rz(-2.3876397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3774776) q[1];
sx q[1];
rz(-0.65122737) q[1];
sx q[1];
rz(2.4098544) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84967126) q[3];
sx q[3];
rz(-2.4115218) q[3];
sx q[3];
rz(1.7770313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98662394) q[2];
sx q[2];
rz(-2.0782317) q[2];
sx q[2];
rz(1.0490136) q[2];
rz(-1.6074041) q[3];
sx q[3];
rz(-1.0020071) q[3];
sx q[3];
rz(-1.775942) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55117115) q[0];
sx q[0];
rz(-0.48749247) q[0];
sx q[0];
rz(-0.76229873) q[0];
rz(0.47362348) q[1];
sx q[1];
rz(-1.381424) q[1];
sx q[1];
rz(2.2332938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.565641) q[0];
sx q[0];
rz(-1.5759908) q[0];
sx q[0];
rz(-2.3495102) q[0];
rz(-pi) q[1];
rz(0.90877675) q[2];
sx q[2];
rz(-2.8637466) q[2];
sx q[2];
rz(1.8709594) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4871074) q[1];
sx q[1];
rz(-1.5562378) q[1];
sx q[1];
rz(1.8076987) q[1];
x q[2];
rz(3.129035) q[3];
sx q[3];
rz(-1.5428484) q[3];
sx q[3];
rz(-1.2716774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0198387) q[2];
sx q[2];
rz(-0.25889954) q[2];
sx q[2];
rz(-0.961595) q[2];
rz(2.614295) q[3];
sx q[3];
rz(-1.3798102) q[3];
sx q[3];
rz(0.56078792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90802646) q[0];
sx q[0];
rz(-2.5416424) q[0];
sx q[0];
rz(-0.34616923) q[0];
rz(-1.2973805) q[1];
sx q[1];
rz(-0.66645122) q[1];
sx q[1];
rz(0.60874879) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52972165) q[0];
sx q[0];
rz(-0.68096113) q[0];
sx q[0];
rz(1.7044071) q[0];
rz(-pi) q[1];
rz(2.4791711) q[2];
sx q[2];
rz(-0.33481827) q[2];
sx q[2];
rz(-1.7551369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.765878) q[1];
sx q[1];
rz(-0.27395136) q[1];
sx q[1];
rz(2.8841444) q[1];
rz(0.56383384) q[3];
sx q[3];
rz(-2.2632416) q[3];
sx q[3];
rz(3.0094224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79494563) q[2];
sx q[2];
rz(-2.5865159) q[2];
sx q[2];
rz(-2.4897599) q[2];
rz(0.73976222) q[3];
sx q[3];
rz(-2.0756105) q[3];
sx q[3];
rz(-2.2447926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5935434) q[0];
sx q[0];
rz(-2.6683922) q[0];
sx q[0];
rz(-2.5439673) q[0];
rz(-1.8400486) q[1];
sx q[1];
rz(-1.5664682) q[1];
sx q[1];
rz(-0.32599932) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6582002) q[0];
sx q[0];
rz(-2.3478386) q[0];
sx q[0];
rz(-1.9824337) q[0];
rz(-pi) q[1];
rz(2.4818899) q[2];
sx q[2];
rz(-2.7355425) q[2];
sx q[2];
rz(-2.7999807) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1797267) q[1];
sx q[1];
rz(-0.71007198) q[1];
sx q[1];
rz(1.1674787) q[1];
x q[2];
rz(1.4526618) q[3];
sx q[3];
rz(-1.6836327) q[3];
sx q[3];
rz(2.0510587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7035383) q[2];
sx q[2];
rz(-1.9376829) q[2];
sx q[2];
rz(2.8013012) q[2];
rz(-1.641364) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(0.59804183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.930645) q[0];
sx q[0];
rz(-1.1177381) q[0];
sx q[0];
rz(-3.0979544) q[0];
rz(-0.71815193) q[1];
sx q[1];
rz(-1.8225881) q[1];
sx q[1];
rz(-1.7125548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315009) q[0];
sx q[0];
rz(-0.25496182) q[0];
sx q[0];
rz(0.73627569) q[0];
x q[1];
rz(-0.67769717) q[2];
sx q[2];
rz(-1.4645828) q[2];
sx q[2];
rz(-0.25523694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2344115) q[1];
sx q[1];
rz(-2.6446807) q[1];
sx q[1];
rz(-0.18053825) q[1];
x q[2];
rz(-1.9472856) q[3];
sx q[3];
rz(-1.2293136) q[3];
sx q[3];
rz(-1.3545413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4622197) q[2];
sx q[2];
rz(-2.2023109) q[2];
sx q[2];
rz(0.77793724) q[2];
rz(-2.098162) q[3];
sx q[3];
rz(-1.5305887) q[3];
sx q[3];
rz(-1.2354318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(3.0935681) q[0];
sx q[0];
rz(-1.0014191) q[0];
sx q[0];
rz(0.77145664) q[0];
rz(1.7780766) q[1];
sx q[1];
rz(-1.172278) q[1];
sx q[1];
rz(-1.5350501) q[1];
rz(-0.58495782) q[2];
sx q[2];
rz(-0.81333209) q[2];
sx q[2];
rz(1.3438136) q[2];
rz(0.89911581) q[3];
sx q[3];
rz(-1.1976783) q[3];
sx q[3];
rz(-1.4730361) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
