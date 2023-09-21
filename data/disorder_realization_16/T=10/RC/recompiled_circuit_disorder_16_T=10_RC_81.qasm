OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(3.9711877) q[0];
sx q[0];
rz(9.2708099) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(0.33831236) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51988039) q[0];
sx q[0];
rz(-1.8147239) q[0];
sx q[0];
rz(-1.2432616) q[0];
rz(-pi) q[1];
rz(-0.36891919) q[2];
sx q[2];
rz(-0.92637617) q[2];
sx q[2];
rz(1.8298139) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0677883) q[1];
sx q[1];
rz(-1.5177625) q[1];
sx q[1];
rz(-2.8552613) q[1];
x q[2];
rz(-1.5987414) q[3];
sx q[3];
rz(-0.51252796) q[3];
sx q[3];
rz(0.41748369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.9677229) q[2];
rz(0.075803444) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(-3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8006111) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(-2.5669572) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(1.2423135) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14318289) q[0];
sx q[0];
rz(-1.3230723) q[0];
sx q[0];
rz(-1.4383016) q[0];
x q[1];
rz(2.9827036) q[2];
sx q[2];
rz(-0.94413589) q[2];
sx q[2];
rz(-1.6275258) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4878792) q[1];
sx q[1];
rz(-0.55760819) q[1];
sx q[1];
rz(-2.8370268) q[1];
rz(-pi) q[2];
rz(-0.55450704) q[3];
sx q[3];
rz(-0.79075659) q[3];
sx q[3];
rz(-2.8046372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80766455) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(-0.58369613) q[2];
rz(-2.5675473) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(-0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7211001) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(0.72845355) q[0];
rz(1.6473673) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(1.0167936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8249614) q[0];
sx q[0];
rz(-1.6065856) q[0];
sx q[0];
rz(-3.0599942) q[0];
rz(-2.3967907) q[2];
sx q[2];
rz(-2.475127) q[2];
sx q[2];
rz(1.9772066) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.46360717) q[1];
sx q[1];
rz(-1.4336168) q[1];
sx q[1];
rz(-0.82171085) q[1];
rz(-3.0924762) q[3];
sx q[3];
rz(-1.2894863) q[3];
sx q[3];
rz(-1.0055055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8016522) q[2];
sx q[2];
rz(-1.5321956) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(-0.63878757) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(-1.1857741) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8124354) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(1.6695492) q[0];
rz(-2.4064348) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(2.8947815) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7128162) q[0];
sx q[0];
rz(-1.4220337) q[0];
sx q[0];
rz(0.87417283) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7814126) q[2];
sx q[2];
rz(-1.1014551) q[2];
sx q[2];
rz(-2.5460555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4775131) q[1];
sx q[1];
rz(-1.2004566) q[1];
sx q[1];
rz(1.5143637) q[1];
x q[2];
rz(-0.28624268) q[3];
sx q[3];
rz(-0.97385588) q[3];
sx q[3];
rz(-0.7569353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(1.5412615) q[2];
rz(0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(2.2533806) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8108869) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(-1.3274308) q[0];
rz(-1.56303) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(0.24838233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7800956) q[0];
sx q[0];
rz(-2.1108315) q[0];
sx q[0];
rz(1.5310775) q[0];
rz(-pi) q[1];
rz(-1.4127172) q[2];
sx q[2];
rz(-1.386596) q[2];
sx q[2];
rz(2.7764111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0548965) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(-2.8413248) q[1];
rz(-pi) q[2];
rz(-1.3316657) q[3];
sx q[3];
rz(-2.2056747) q[3];
sx q[3];
rz(2.2374416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7053232) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(2.3441337) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(1.3457993) q[0];
rz(1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(-3.016901) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7170982) q[0];
sx q[0];
rz(-1.4846804) q[0];
sx q[0];
rz(-2.9647102) q[0];
rz(-pi) q[1];
rz(-2.7331946) q[2];
sx q[2];
rz(-0.64986594) q[2];
sx q[2];
rz(1.0330531) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71619892) q[1];
sx q[1];
rz(-0.73912207) q[1];
sx q[1];
rz(-3.0576586) q[1];
rz(-2.3966339) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(1.484364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51320118) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(-2.0882873) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(-2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597647) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(2.0986309) q[0];
rz(0.46328059) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(-1.0707062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28744222) q[0];
sx q[0];
rz(-1.5713912) q[0];
sx q[0];
rz(-2.6575412) q[0];
x q[1];
rz(-1.7037017) q[2];
sx q[2];
rz(-1.9153567) q[2];
sx q[2];
rz(0.29495707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6003905) q[1];
sx q[1];
rz(-1.6704847) q[1];
sx q[1];
rz(0.77460918) q[1];
x q[2];
rz(-1.5504863) q[3];
sx q[3];
rz(-1.9029402) q[3];
sx q[3];
rz(0.16196812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36859194) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-0.32361844) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(-1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535646) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(-1.2063684) q[0];
rz(1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(2.1941197) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5652126) q[0];
sx q[0];
rz(-1.0351666) q[0];
sx q[0];
rz(-3.0239848) q[0];
x q[1];
rz(-2.410789) q[2];
sx q[2];
rz(-0.23554221) q[2];
sx q[2];
rz(-1.8159602) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35093388) q[1];
sx q[1];
rz(-2.5194063) q[1];
sx q[1];
rz(1.15637) q[1];
x q[2];
rz(0.8074699) q[3];
sx q[3];
rz(-1.315457) q[3];
sx q[3];
rz(1.0367928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5802713) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(0.46978152) q[2];
rz(-1.3011159) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(-2.8619213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(1.8126194) q[0];
rz(0.7912311) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(2.9387617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79265362) q[0];
sx q[0];
rz(-1.5968423) q[0];
sx q[0];
rz(0.039020122) q[0];
rz(0.22325309) q[2];
sx q[2];
rz(-1.7749783) q[2];
sx q[2];
rz(-2.5399361) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.030414) q[1];
sx q[1];
rz(-1.5215538) q[1];
sx q[1];
rz(-1.6715675) q[1];
rz(-pi) q[2];
rz(0.46189724) q[3];
sx q[3];
rz(-1.4326722) q[3];
sx q[3];
rz(-1.671333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(-2.1450796) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-2.3341808) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0614232) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(-2.9598575) q[0];
rz(3.0985447) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(-2.8607686) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.735606) q[0];
sx q[0];
rz(-1.1445023) q[0];
sx q[0];
rz(2.6398525) q[0];
rz(-2.4091987) q[2];
sx q[2];
rz(-2.8576982) q[2];
sx q[2];
rz(2.4925799) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9210789) q[1];
sx q[1];
rz(-1.5346569) q[1];
sx q[1];
rz(-1.6284579) q[1];
x q[2];
rz(1.9032352) q[3];
sx q[3];
rz(-1.7566163) q[3];
sx q[3];
rz(-0.18248724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8250371) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(0.62310702) q[2];
rz(2.1394219) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.4939209) q[0];
sx q[0];
rz(-1.5734084) q[0];
sx q[0];
rz(-1.5403803) q[0];
rz(-0.87396809) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-0.95332425) q[2];
sx q[2];
rz(-0.8669903) q[2];
sx q[2];
rz(0.16624761) q[2];
rz(-1.9813886) q[3];
sx q[3];
rz(-1.900233) q[3];
sx q[3];
rz(-1.486447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];