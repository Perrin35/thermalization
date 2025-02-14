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
rz(1.7841568) q[0];
sx q[0];
rz(-0.62499243) q[0];
sx q[0];
rz(1.5224737) q[0];
rz(-1.9982665) q[1];
sx q[1];
rz(-0.63653094) q[1];
sx q[1];
rz(-1.7826537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1899043) q[0];
sx q[0];
rz(-1.0190599) q[0];
sx q[0];
rz(-0.34073982) q[0];
x q[1];
rz(-1.9404561) q[2];
sx q[2];
rz(-1.2944752) q[2];
sx q[2];
rz(-1.2075961) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93329079) q[1];
sx q[1];
rz(-1.8027026) q[1];
sx q[1];
rz(-1.8342706) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10349689) q[3];
sx q[3];
rz(-2.0465133) q[3];
sx q[3];
rz(-0.58799839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6978567) q[2];
sx q[2];
rz(-2.3190658) q[2];
sx q[2];
rz(-0.47625429) q[2];
rz(3.0016628) q[3];
sx q[3];
rz(-1.5014476) q[3];
sx q[3];
rz(3.068315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.78854617) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(2.7781558) q[0];
rz(0.67047554) q[1];
sx q[1];
rz(-1.8164219) q[1];
sx q[1];
rz(1.0796116) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19980656) q[0];
sx q[0];
rz(-2.5528756) q[0];
sx q[0];
rz(-1.1325645) q[0];
rz(-pi) q[1];
x q[1];
rz(3.011376) q[2];
sx q[2];
rz(-1.0417582) q[2];
sx q[2];
rz(1.0769653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.590821) q[1];
sx q[1];
rz(-1.5268004) q[1];
sx q[1];
rz(-1.8122871) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5609301) q[3];
sx q[3];
rz(-1.7911909) q[3];
sx q[3];
rz(3.0739258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67538992) q[2];
sx q[2];
rz(-0.73489302) q[2];
sx q[2];
rz(-2.7030763) q[2];
rz(-2.130326) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(-1.6605759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48280516) q[0];
sx q[0];
rz(-2.6130982) q[0];
sx q[0];
rz(2.3129789) q[0];
rz(1.2476745) q[1];
sx q[1];
rz(-1.1670466) q[1];
sx q[1];
rz(0.25965986) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22406507) q[0];
sx q[0];
rz(-0.7389141) q[0];
sx q[0];
rz(-1.1631484) q[0];
rz(-0.38085085) q[2];
sx q[2];
rz(-1.9162188) q[2];
sx q[2];
rz(2.7598515) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6344172) q[1];
sx q[1];
rz(-0.89752642) q[1];
sx q[1];
rz(-2.5310764) q[1];
x q[2];
rz(-2.8494495) q[3];
sx q[3];
rz(-2.202987) q[3];
sx q[3];
rz(1.2911673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9344249) q[2];
sx q[2];
rz(-1.577689) q[2];
sx q[2];
rz(-0.820532) q[2];
rz(0.82646838) q[3];
sx q[3];
rz(-0.99244899) q[3];
sx q[3];
rz(-1.8291399) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25927037) q[0];
sx q[0];
rz(-2.3311908) q[0];
sx q[0];
rz(0.93233863) q[0];
rz(1.3177634) q[1];
sx q[1];
rz(-1.1612786) q[1];
sx q[1];
rz(3.0188149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6500191) q[0];
sx q[0];
rz(-0.54560018) q[0];
sx q[0];
rz(0.16262098) q[0];
x q[1];
rz(0.30722801) q[2];
sx q[2];
rz(-0.84103742) q[2];
sx q[2];
rz(0.73452362) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.40298324) q[1];
sx q[1];
rz(-1.5697641) q[1];
sx q[1];
rz(0.82768517) q[1];
rz(-pi) q[2];
rz(2.6969321) q[3];
sx q[3];
rz(-0.61189368) q[3];
sx q[3];
rz(-1.1763619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54470283) q[2];
sx q[2];
rz(-2.2201316) q[2];
sx q[2];
rz(0.60453647) q[2];
rz(0.70945159) q[3];
sx q[3];
rz(-2.3716898) q[3];
sx q[3];
rz(-0.82367045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17968793) q[0];
sx q[0];
rz(-3.0178495) q[0];
sx q[0];
rz(-1.3667579) q[0];
rz(2.1429515) q[1];
sx q[1];
rz(-2.3256681) q[1];
sx q[1];
rz(0.5407812) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3872469) q[0];
sx q[0];
rz(-1.7779113) q[0];
sx q[0];
rz(-1.4282754) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2774463) q[2];
sx q[2];
rz(-0.96071488) q[2];
sx q[2];
rz(0.44877258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.97637002) q[1];
sx q[1];
rz(-1.762855) q[1];
sx q[1];
rz(-2.8403175) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56638797) q[3];
sx q[3];
rz(-1.4175804) q[3];
sx q[3];
rz(-0.18463102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46405408) q[2];
sx q[2];
rz(-2.3276734) q[2];
sx q[2];
rz(2.3133254) q[2];
rz(-1.6895435) q[3];
sx q[3];
rz(-1.6477081) q[3];
sx q[3];
rz(0.90788466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.117711) q[0];
sx q[0];
rz(-1.9458867) q[0];
sx q[0];
rz(-1.1718132) q[0];
rz(-2.4614629) q[1];
sx q[1];
rz(-1.9155733) q[1];
sx q[1];
rz(0.5562869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2883462) q[0];
sx q[0];
rz(-0.46015209) q[0];
sx q[0];
rz(0.96921885) q[0];
x q[1];
rz(0.0027782253) q[2];
sx q[2];
rz(-1.7173212) q[2];
sx q[2];
rz(-0.79039449) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57842931) q[1];
sx q[1];
rz(-0.99908864) q[1];
sx q[1];
rz(1.494638) q[1];
rz(-pi) q[2];
rz(-0.63594922) q[3];
sx q[3];
rz(-1.9665641) q[3];
sx q[3];
rz(1.1102939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3351626) q[2];
sx q[2];
rz(-1.282629) q[2];
sx q[2];
rz(2.9365017) q[2];
rz(2.3388376) q[3];
sx q[3];
rz(-1.1057248) q[3];
sx q[3];
rz(2.5870489) q[3];
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
rz(-pi) q[0];
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
rz(-0.71451181) q[0];
sx q[0];
rz(-1.0733805) q[0];
sx q[0];
rz(0.22739534) q[0];
rz(2.3511476) q[1];
sx q[1];
rz(-0.77722725) q[1];
sx q[1];
rz(-2.3048185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8604383) q[0];
sx q[0];
rz(-1.3520762) q[0];
sx q[0];
rz(2.6118096) q[0];
rz(-1.2143986) q[2];
sx q[2];
rz(-2.1292392) q[2];
sx q[2];
rz(1.7624378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62722271) q[1];
sx q[1];
rz(-0.52241507) q[1];
sx q[1];
rz(0.66607968) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3804139) q[3];
sx q[3];
rz(-1.8463148) q[3];
sx q[3];
rz(2.8305284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5922015) q[2];
sx q[2];
rz(-1.6624007) q[2];
sx q[2];
rz(0.57360348) q[2];
rz(-2.8163689) q[3];
sx q[3];
rz(-0.89541382) q[3];
sx q[3];
rz(-0.4044683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10930546) q[0];
sx q[0];
rz(-2.9635297) q[0];
sx q[0];
rz(0.67894116) q[0];
rz(-3.0067054) q[1];
sx q[1];
rz(-2.0811681) q[1];
sx q[1];
rz(-1.7178242) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7622691) q[0];
sx q[0];
rz(-1.0417308) q[0];
sx q[0];
rz(2.6075415) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5871954) q[2];
sx q[2];
rz(-2.1986678) q[2];
sx q[2];
rz(0.45931057) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.4383761) q[1];
sx q[1];
rz(-1.5408311) q[1];
sx q[1];
rz(2.334146) q[1];
x q[2];
rz(2.0349488) q[3];
sx q[3];
rz(-1.8309793) q[3];
sx q[3];
rz(1.2342208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5474825) q[2];
sx q[2];
rz(-0.91071931) q[2];
sx q[2];
rz(2.9098848) q[2];
rz(-1.0101275) q[3];
sx q[3];
rz(-1.908327) q[3];
sx q[3];
rz(-2.2505984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4657087) q[0];
sx q[0];
rz(-0.75513419) q[0];
sx q[0];
rz(-0.51374197) q[0];
rz(1.7087917) q[1];
sx q[1];
rz(-1.2756196) q[1];
sx q[1];
rz(-1.6339711) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20376539) q[0];
sx q[0];
rz(-1.7047593) q[0];
sx q[0];
rz(-1.6197657) q[0];
rz(2.5812923) q[2];
sx q[2];
rz(-2.8447084) q[2];
sx q[2];
rz(0.96499824) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0158314) q[1];
sx q[1];
rz(-1.4471701) q[1];
sx q[1];
rz(3.0112253) q[1];
rz(-1.8225372) q[3];
sx q[3];
rz(-1.5949524) q[3];
sx q[3];
rz(2.5509953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.12019176) q[2];
sx q[2];
rz(-2.2970565) q[2];
sx q[2];
rz(2.8821442) q[2];
rz(-1.1546968) q[3];
sx q[3];
rz(-0.76609937) q[3];
sx q[3];
rz(1.2999138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65086377) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(1.65253) q[0];
rz(-2.5015855) q[1];
sx q[1];
rz(-2.044544) q[1];
sx q[1];
rz(0.99001137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31371597) q[0];
sx q[0];
rz(-1.6589239) q[0];
sx q[0];
rz(-0.14493305) q[0];
rz(-2.0406538) q[2];
sx q[2];
rz(-0.85950101) q[2];
sx q[2];
rz(2.2180722) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.99553019) q[1];
sx q[1];
rz(-1.4096673) q[1];
sx q[1];
rz(-1.9934158) q[1];
rz(1.680671) q[3];
sx q[3];
rz(-0.2239033) q[3];
sx q[3];
rz(-2.2810625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2085569) q[2];
sx q[2];
rz(-1.8427589) q[2];
sx q[2];
rz(-0.2253069) q[2];
rz(0.92132583) q[3];
sx q[3];
rz(-2.7511629) q[3];
sx q[3];
rz(3.0909753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13023547) q[0];
sx q[0];
rz(-0.36778944) q[0];
sx q[0];
rz(0.854048) q[0];
rz(1.8232952) q[1];
sx q[1];
rz(-1.6164936) q[1];
sx q[1];
rz(2.2093538) q[1];
rz(0.097250289) q[2];
sx q[2];
rz(-1.8293867) q[2];
sx q[2];
rz(0.23501227) q[2];
rz(0.18410889) q[3];
sx q[3];
rz(-2.1118023) q[3];
sx q[3];
rz(2.8258256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
