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
rz(0.99474466) q[0];
sx q[0];
rz(3.7288546) q[0];
sx q[0];
rz(10.050339) q[0];
rz(0.62043959) q[1];
sx q[1];
rz(-2.151139) q[1];
sx q[1];
rz(0.78392309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44255689) q[0];
sx q[0];
rz(-0.43414206) q[0];
sx q[0];
rz(0.0084369466) q[0];
x q[1];
rz(3.002499) q[2];
sx q[2];
rz(-1.6533829) q[2];
sx q[2];
rz(-0.54143426) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3202438) q[1];
sx q[1];
rz(-0.87961266) q[1];
sx q[1];
rz(-2.3752179) q[1];
rz(2.3638681) q[3];
sx q[3];
rz(-1.6334264) q[3];
sx q[3];
rz(0.81111139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1938532) q[2];
sx q[2];
rz(-1.0449907) q[2];
sx q[2];
rz(-0.6673153) q[2];
rz(0.31467485) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(-0.92710322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84745234) q[0];
sx q[0];
rz(-0.84646928) q[0];
sx q[0];
rz(0.2997998) q[0];
rz(1.0046129) q[1];
sx q[1];
rz(-0.71669465) q[1];
sx q[1];
rz(-0.85003781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4097828) q[0];
sx q[0];
rz(-0.96854051) q[0];
sx q[0];
rz(0.11645634) q[0];
x q[1];
rz(0.17053594) q[2];
sx q[2];
rz(-1.6493622) q[2];
sx q[2];
rz(-2.6002392) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0850573) q[1];
sx q[1];
rz(-1.3621482) q[1];
sx q[1];
rz(1.2056808) q[1];
rz(1.5845234) q[3];
sx q[3];
rz(-1.1662332) q[3];
sx q[3];
rz(-1.3409529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5726418) q[2];
sx q[2];
rz(-2.6093542) q[2];
sx q[2];
rz(2.4786095) q[2];
rz(0.7524544) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74333423) q[0];
sx q[0];
rz(-2.7404116) q[0];
sx q[0];
rz(-2.4617526) q[0];
rz(2.4196449) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(-0.78071761) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762622) q[0];
sx q[0];
rz(-1.1701382) q[0];
sx q[0];
rz(-0.0082992502) q[0];
rz(-pi) q[1];
rz(-1.5268097) q[2];
sx q[2];
rz(-3.0553748) q[2];
sx q[2];
rz(1.9808863) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5131849) q[1];
sx q[1];
rz(-2.3806913) q[1];
sx q[1];
rz(0.55879467) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3296919) q[3];
sx q[3];
rz(-0.74987292) q[3];
sx q[3];
rz(-2.1151224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9146933) q[2];
sx q[2];
rz(-1.2279899) q[2];
sx q[2];
rz(-2.4719888) q[2];
rz(-2.2412444) q[3];
sx q[3];
rz(-2.1293631) q[3];
sx q[3];
rz(-0.77696925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8825619) q[0];
sx q[0];
rz(-2.7029523) q[0];
sx q[0];
rz(-0.90748179) q[0];
rz(2.5471845) q[1];
sx q[1];
rz(-0.92955697) q[1];
sx q[1];
rz(-2.0118735) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0967207) q[0];
sx q[0];
rz(-1.7611928) q[0];
sx q[0];
rz(-2.2119766) q[0];
rz(-pi) q[1];
rz(1.8215034) q[2];
sx q[2];
rz(-0.54566979) q[2];
sx q[2];
rz(-2.9755993) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10953021) q[1];
sx q[1];
rz(-2.1550165) q[1];
sx q[1];
rz(-1.24448) q[1];
rz(-0.28092556) q[3];
sx q[3];
rz(-1.2884029) q[3];
sx q[3];
rz(0.97746935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0323459) q[2];
sx q[2];
rz(-1.7416019) q[2];
sx q[2];
rz(-2.4692811) q[2];
rz(-0.0191056) q[3];
sx q[3];
rz(-0.34135434) q[3];
sx q[3];
rz(-0.13489558) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30447176) q[0];
sx q[0];
rz(-2.3470375) q[0];
sx q[0];
rz(-0.16803148) q[0];
rz(-1.9145603) q[1];
sx q[1];
rz(-1.2201759) q[1];
sx q[1];
rz(2.8438445) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79556161) q[0];
sx q[0];
rz(-2.3712284) q[0];
sx q[0];
rz(2.3576231) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3565294) q[2];
sx q[2];
rz(-0.74399555) q[2];
sx q[2];
rz(0.97563484) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9492901) q[1];
sx q[1];
rz(-1.5459035) q[1];
sx q[1];
rz(-0.20583238) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7295804) q[3];
sx q[3];
rz(-1.571469) q[3];
sx q[3];
rz(-2.1070045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1905404) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(0.61868787) q[2];
rz(-0.80858532) q[3];
sx q[3];
rz(-2.6638668) q[3];
sx q[3];
rz(1.3396858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65281463) q[0];
sx q[0];
rz(-1.6365016) q[0];
sx q[0];
rz(0.46522796) q[0];
rz(0.48509994) q[1];
sx q[1];
rz(-0.41050375) q[1];
sx q[1];
rz(3.0533275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34644914) q[0];
sx q[0];
rz(-0.69899054) q[0];
sx q[0];
rz(-0.82620746) q[0];
rz(1.5851998) q[2];
sx q[2];
rz(-0.72897899) q[2];
sx q[2];
rz(0.62504238) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3151206) q[1];
sx q[1];
rz(-1.406028) q[1];
sx q[1];
rz(-0.75854782) q[1];
rz(-1.4824105) q[3];
sx q[3];
rz(-1.2840349) q[3];
sx q[3];
rz(0.6500611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39468592) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(0.49760231) q[2];
rz(-0.48007128) q[3];
sx q[3];
rz(-2.2205455) q[3];
sx q[3];
rz(-3.0729496) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29480544) q[0];
sx q[0];
rz(-0.46600309) q[0];
sx q[0];
rz(1.4303327) q[0];
rz(-1.826674) q[1];
sx q[1];
rz(-0.39350915) q[1];
sx q[1];
rz(1.5865145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85077943) q[0];
sx q[0];
rz(-1.4516136) q[0];
sx q[0];
rz(1.9656154) q[0];
rz(-pi) q[1];
rz(2.2847963) q[2];
sx q[2];
rz(-1.1044377) q[2];
sx q[2];
rz(-3.015632) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79826971) q[1];
sx q[1];
rz(-2.1499691) q[1];
sx q[1];
rz(-0.12076853) q[1];
x q[2];
rz(-0.93933834) q[3];
sx q[3];
rz(-1.5087481) q[3];
sx q[3];
rz(-2.5417059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.84307182) q[2];
sx q[2];
rz(-2.4303747) q[2];
sx q[2];
rz(2.2141875) q[2];
rz(-1.983042) q[3];
sx q[3];
rz(-2.4535593) q[3];
sx q[3];
rz(0.096175171) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226828) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(-0.476015) q[0];
rz(2.1503275) q[1];
sx q[1];
rz(-2.3920993) q[1];
sx q[1];
rz(0.083866619) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96362417) q[0];
sx q[0];
rz(-1.2126847) q[0];
sx q[0];
rz(-2.0923418) q[0];
rz(-pi) q[1];
rz(1.2071361) q[2];
sx q[2];
rz(-2.0061099) q[2];
sx q[2];
rz(2.2029623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8955651) q[1];
sx q[1];
rz(-1.9894587) q[1];
sx q[1];
rz(3.0484867) q[1];
rz(-3.0552577) q[3];
sx q[3];
rz(-1.9342967) q[3];
sx q[3];
rz(0.32237651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91741651) q[2];
sx q[2];
rz(-0.26599628) q[2];
sx q[2];
rz(-2.4931397) q[2];
rz(0.90483856) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(0.31074935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7998578) q[0];
sx q[0];
rz(-2.4315727) q[0];
sx q[0];
rz(-0.56890666) q[0];
rz(2.3078602) q[1];
sx q[1];
rz(-1.2833475) q[1];
sx q[1];
rz(-2.1368829) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6093401) q[0];
sx q[0];
rz(-0.94826094) q[0];
sx q[0];
rz(-1.3329766) q[0];
x q[1];
rz(1.0332074) q[2];
sx q[2];
rz(-1.7694607) q[2];
sx q[2];
rz(-0.63420701) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9819599) q[1];
sx q[1];
rz(-1.4985604) q[1];
sx q[1];
rz(0.87149974) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58752693) q[3];
sx q[3];
rz(-1.5024795) q[3];
sx q[3];
rz(-0.63623601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4363165) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(0.3479859) q[2];
rz(-0.82585382) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(-0.60514778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0843049) q[0];
sx q[0];
rz(-1.0658406) q[0];
sx q[0];
rz(-2.9344946) q[0];
rz(-2.6026978) q[1];
sx q[1];
rz(-2.5205044) q[1];
sx q[1];
rz(-2.5394687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3054863) q[0];
sx q[0];
rz(-2.0417111) q[0];
sx q[0];
rz(-2.6013646) q[0];
rz(-pi) q[1];
rz(-0.57453491) q[2];
sx q[2];
rz(-2.4705187) q[2];
sx q[2];
rz(0.35071638) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.94127266) q[1];
sx q[1];
rz(-1.4956022) q[1];
sx q[1];
rz(1.978051) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8982417) q[3];
sx q[3];
rz(-0.85059887) q[3];
sx q[3];
rz(0.44046381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78486097) q[2];
sx q[2];
rz(-0.95180231) q[2];
sx q[2];
rz(-1.7968563) q[2];
rz(0.52062672) q[3];
sx q[3];
rz(-0.27025637) q[3];
sx q[3];
rz(-2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44705924) q[0];
sx q[0];
rz(-2.4166528) q[0];
sx q[0];
rz(1.7550533) q[0];
rz(0.78320349) q[1];
sx q[1];
rz(-1.422311) q[1];
sx q[1];
rz(1.5625988) q[1];
rz(-1.3010018) q[2];
sx q[2];
rz(-1.9023583) q[2];
sx q[2];
rz(-2.6826774) q[2];
rz(-0.047719638) q[3];
sx q[3];
rz(-0.24038355) q[3];
sx q[3];
rz(3.0265831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
