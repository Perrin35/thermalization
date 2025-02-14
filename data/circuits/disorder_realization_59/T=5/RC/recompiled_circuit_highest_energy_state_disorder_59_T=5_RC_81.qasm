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
rz(0.94035971) q[0];
sx q[0];
rz(-2.0191963) q[0];
sx q[0];
rz(2.9474592) q[0];
rz(1.0951207) q[1];
sx q[1];
rz(-1.6835338) q[1];
sx q[1];
rz(-0.74498743) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783837) q[0];
sx q[0];
rz(-1.4683405) q[0];
sx q[0];
rz(-1.6215023) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81211798) q[2];
sx q[2];
rz(-1.0924763) q[2];
sx q[2];
rz(1.1253192) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2304334) q[1];
sx q[1];
rz(-1.3285561) q[1];
sx q[1];
rz(-1.469286) q[1];
rz(-0.93985112) q[3];
sx q[3];
rz(-0.54169023) q[3];
sx q[3];
rz(1.9909137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3377043) q[2];
sx q[2];
rz(-1.3362198) q[2];
sx q[2];
rz(1.3275006) q[2];
rz(-1.368807) q[3];
sx q[3];
rz(-0.33392206) q[3];
sx q[3];
rz(-2.9201065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9889481) q[0];
sx q[0];
rz(-1.7247609) q[0];
sx q[0];
rz(-0.055135559) q[0];
rz(-2.4368743) q[1];
sx q[1];
rz(-1.5635419) q[1];
sx q[1];
rz(1.0447186) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9023414) q[0];
sx q[0];
rz(-1.2980868) q[0];
sx q[0];
rz(-0.38614778) q[0];
rz(2.3364445) q[2];
sx q[2];
rz(-2.3986446) q[2];
sx q[2];
rz(0.93490237) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5767908) q[1];
sx q[1];
rz(-1.163244) q[1];
sx q[1];
rz(2.3840586) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0468076) q[3];
sx q[3];
rz(-1.6322502) q[3];
sx q[3];
rz(-0.8549785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3220871) q[2];
sx q[2];
rz(-0.6531859) q[2];
sx q[2];
rz(1.8196003) q[2];
rz(-0.78754887) q[3];
sx q[3];
rz(-1.4044263) q[3];
sx q[3];
rz(-1.0549841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2128485) q[0];
sx q[0];
rz(-0.66900122) q[0];
sx q[0];
rz(-2.4231732) q[0];
rz(-2.8058715) q[1];
sx q[1];
rz(-1.4664783) q[1];
sx q[1];
rz(2.1885923) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9586784) q[0];
sx q[0];
rz(-1.0299068) q[0];
sx q[0];
rz(3.0989585) q[0];
rz(-0.14710958) q[2];
sx q[2];
rz(-2.0858371) q[2];
sx q[2];
rz(1.7619606) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6086207) q[1];
sx q[1];
rz(-0.91580456) q[1];
sx q[1];
rz(2.6650732) q[1];
x q[2];
rz(0.6742874) q[3];
sx q[3];
rz(-2.4113048) q[3];
sx q[3];
rz(-1.3004608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0042051729) q[2];
sx q[2];
rz(-1.1342099) q[2];
sx q[2];
rz(0.65286621) q[2];
rz(2.3434434) q[3];
sx q[3];
rz(-1.8316725) q[3];
sx q[3];
rz(0.64364141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2586486) q[0];
sx q[0];
rz(-0.84027165) q[0];
sx q[0];
rz(0.81664455) q[0];
rz(-0.69149292) q[1];
sx q[1];
rz(-1.34812) q[1];
sx q[1];
rz(0.74849558) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8832421) q[0];
sx q[0];
rz(-1.878698) q[0];
sx q[0];
rz(1.7791722) q[0];
rz(-pi) q[1];
rz(0.091684999) q[2];
sx q[2];
rz(-1.3576512) q[2];
sx q[2];
rz(-0.64848778) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2556127) q[1];
sx q[1];
rz(-1.1062262) q[1];
sx q[1];
rz(0.8452615) q[1];
rz(-0.17552142) q[3];
sx q[3];
rz(-1.3427882) q[3];
sx q[3];
rz(1.0358178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3889435) q[2];
sx q[2];
rz(-0.77072531) q[2];
sx q[2];
rz(0.79461092) q[2];
rz(1.4258344) q[3];
sx q[3];
rz(-1.040193) q[3];
sx q[3];
rz(0.37253255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36885095) q[0];
sx q[0];
rz(-2.0617101) q[0];
sx q[0];
rz(-2.6439731) q[0];
rz(2.3720062) q[1];
sx q[1];
rz(-1.9895357) q[1];
sx q[1];
rz(2.41113) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896547) q[0];
sx q[0];
rz(-1.7853072) q[0];
sx q[0];
rz(0.69846054) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1986352) q[2];
sx q[2];
rz(-1.3636317) q[2];
sx q[2];
rz(3.0414273) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.5268908) q[1];
sx q[1];
rz(-2.3040471) q[1];
sx q[1];
rz(2.8575543) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3114871) q[3];
sx q[3];
rz(-2.5200994) q[3];
sx q[3];
rz(1.9754387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2490425) q[2];
sx q[2];
rz(-1.5183134) q[2];
sx q[2];
rz(2.6882233) q[2];
rz(1.310965) q[3];
sx q[3];
rz(-0.1736719) q[3];
sx q[3];
rz(3.0176676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740358) q[0];
sx q[0];
rz(-2.0581764) q[0];
sx q[0];
rz(1.3483082) q[0];
rz(2.0454316) q[1];
sx q[1];
rz(-1.4827012) q[1];
sx q[1];
rz(-0.62166628) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28866119) q[0];
sx q[0];
rz(-0.71446361) q[0];
sx q[0];
rz(3.1174401) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6699583) q[2];
sx q[2];
rz(-1.1767724) q[2];
sx q[2];
rz(0.1876012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7335947) q[1];
sx q[1];
rz(-0.80469614) q[1];
sx q[1];
rz(0.3478653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4824941) q[3];
sx q[3];
rz(-2.3602398) q[3];
sx q[3];
rz(0.78418651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93214503) q[2];
sx q[2];
rz(-2.3422082) q[2];
sx q[2];
rz(1.4383593) q[2];
rz(0.98073331) q[3];
sx q[3];
rz(-1.8756198) q[3];
sx q[3];
rz(-1.929662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.257523) q[0];
sx q[0];
rz(-1.3781837) q[0];
sx q[0];
rz(2.8164465) q[0];
rz(2.6349321) q[1];
sx q[1];
rz(-1.0136565) q[1];
sx q[1];
rz(-1.8258757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36111212) q[0];
sx q[0];
rz(-0.69421116) q[0];
sx q[0];
rz(-1.3785465) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99847128) q[2];
sx q[2];
rz(-1.9461759) q[2];
sx q[2];
rz(1.4294415) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.084242227) q[1];
sx q[1];
rz(-2.2154543) q[1];
sx q[1];
rz(1.0727543) q[1];
rz(-pi) q[2];
rz(-2.2030284) q[3];
sx q[3];
rz(-1.5364416) q[3];
sx q[3];
rz(-1.1704579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.331984) q[2];
sx q[2];
rz(-2.1029682) q[2];
sx q[2];
rz(-1.6161551) q[2];
rz(-1.7841313) q[3];
sx q[3];
rz(-1.6922035) q[3];
sx q[3];
rz(-1.6239032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9759393) q[0];
sx q[0];
rz(-1.140056) q[0];
sx q[0];
rz(-2.4923988) q[0];
rz(0.80750418) q[1];
sx q[1];
rz(-2.0048001) q[1];
sx q[1];
rz(1.3884707) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20304414) q[0];
sx q[0];
rz(-1.6106092) q[0];
sx q[0];
rz(0.96588246) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0031707) q[2];
sx q[2];
rz(-1.8571203) q[2];
sx q[2];
rz(0.22655205) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5877046) q[1];
sx q[1];
rz(-2.9575051) q[1];
sx q[1];
rz(-1.192993) q[1];
x q[2];
rz(0.3921359) q[3];
sx q[3];
rz(-1.3754505) q[3];
sx q[3];
rz(2.3980262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93747741) q[2];
sx q[2];
rz(-2.9948339) q[2];
sx q[2];
rz(-2.333763) q[2];
rz(-2.4701123) q[3];
sx q[3];
rz(-1.5059794) q[3];
sx q[3];
rz(1.5427264) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2272187) q[0];
sx q[0];
rz(-2.6052573) q[0];
sx q[0];
rz(2.6887023) q[0];
rz(2.8462786) q[1];
sx q[1];
rz(-1.2295877) q[1];
sx q[1];
rz(-1.3989075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5495396) q[0];
sx q[0];
rz(-0.73707923) q[0];
sx q[0];
rz(1.6960742) q[0];
rz(1.2981425) q[2];
sx q[2];
rz(-2.6116707) q[2];
sx q[2];
rz(0.31604813) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4860515) q[1];
sx q[1];
rz(-0.60645804) q[1];
sx q[1];
rz(0.70041366) q[1];
rz(-pi) q[2];
rz(-1.3736428) q[3];
sx q[3];
rz(-1.592336) q[3];
sx q[3];
rz(-2.7113999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7704775) q[2];
sx q[2];
rz(-0.8728084) q[2];
sx q[2];
rz(-2.5650043) q[2];
rz(-2.9234486) q[3];
sx q[3];
rz(-2.348867) q[3];
sx q[3];
rz(1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169154) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(0.088454811) q[0];
rz(-0.5312008) q[1];
sx q[1];
rz(-0.4041268) q[1];
sx q[1];
rz(-1.6204576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3192465) q[0];
sx q[0];
rz(-0.40238956) q[0];
sx q[0];
rz(-1.7850384) q[0];
rz(-1.0934866) q[2];
sx q[2];
rz(-1.6132406) q[2];
sx q[2];
rz(2.1316656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0853569) q[1];
sx q[1];
rz(-1.7952982) q[1];
sx q[1];
rz(1.4957929) q[1];
rz(0.59595396) q[3];
sx q[3];
rz(-1.9203382) q[3];
sx q[3];
rz(-0.28448018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8497808) q[2];
sx q[2];
rz(-1.3498638) q[2];
sx q[2];
rz(-1.6144217) q[2];
rz(-1.8725066) q[3];
sx q[3];
rz(-0.98617712) q[3];
sx q[3];
rz(1.0482949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0328746) q[0];
sx q[0];
rz(-2.0357108) q[0];
sx q[0];
rz(2.3172814) q[0];
rz(0.4737919) q[1];
sx q[1];
rz(-1.5693249) q[1];
sx q[1];
rz(1.5798461) q[1];
rz(-0.5886336) q[2];
sx q[2];
rz(-2.4559173) q[2];
sx q[2];
rz(0.46951175) q[2];
rz(0.99226034) q[3];
sx q[3];
rz(-1.0836061) q[3];
sx q[3];
rz(1.5988812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
