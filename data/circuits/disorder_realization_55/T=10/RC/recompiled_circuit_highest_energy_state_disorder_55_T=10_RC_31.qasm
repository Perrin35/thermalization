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
rz(3.6291549) q[0];
sx q[0];
rz(3.5222375) q[0];
sx q[0];
rz(9.335523) q[0];
rz(2.6788977) q[1];
sx q[1];
rz(-0.2806288) q[1];
sx q[1];
rz(1.3624924) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87088764) q[0];
sx q[0];
rz(-1.9438367) q[0];
sx q[0];
rz(-1.5482658) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5342185) q[2];
sx q[2];
rz(-1.4223411) q[2];
sx q[2];
rz(-2.0853993) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57682788) q[1];
sx q[1];
rz(-0.62678465) q[1];
sx q[1];
rz(-0.38615055) q[1];
rz(-2.9611601) q[3];
sx q[3];
rz(-2.6789634) q[3];
sx q[3];
rz(-2.608046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6224299) q[2];
sx q[2];
rz(-2.3795655) q[2];
sx q[2];
rz(0.92987522) q[2];
rz(-0.30396384) q[3];
sx q[3];
rz(-1.8094939) q[3];
sx q[3];
rz(2.9728594) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836483) q[0];
sx q[0];
rz(-1.5771663) q[0];
sx q[0];
rz(-1.8899348) q[0];
rz(2.9257863) q[1];
sx q[1];
rz(-2.1025175) q[1];
sx q[1];
rz(-3.0024517) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3045791) q[0];
sx q[0];
rz(-2.4615335) q[0];
sx q[0];
rz(1.3152298) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6670447) q[2];
sx q[2];
rz(-0.54891551) q[2];
sx q[2];
rz(2.7868556) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5522668) q[1];
sx q[1];
rz(-1.8084452) q[1];
sx q[1];
rz(2.4219805) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4226133) q[3];
sx q[3];
rz(-0.96299473) q[3];
sx q[3];
rz(1.3083576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6824049) q[2];
sx q[2];
rz(-0.46899691) q[2];
sx q[2];
rz(-2.0501308) q[2];
rz(1.5952716) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(2.6923164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062155398) q[0];
sx q[0];
rz(-1.9154444) q[0];
sx q[0];
rz(-3.1297744) q[0];
rz(-0.91105175) q[1];
sx q[1];
rz(-1.751519) q[1];
sx q[1];
rz(-1.8131088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8810711) q[0];
sx q[0];
rz(-1.0539454) q[0];
sx q[0];
rz(-1.3840186) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19085796) q[2];
sx q[2];
rz(-1.4578739) q[2];
sx q[2];
rz(2.0115122) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.40483002) q[1];
sx q[1];
rz(-1.9556082) q[1];
sx q[1];
rz(-0.75508187) q[1];
rz(-0.99081466) q[3];
sx q[3];
rz(-3.0085251) q[3];
sx q[3];
rz(-1.7855438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38982424) q[2];
sx q[2];
rz(-1.5247034) q[2];
sx q[2];
rz(2.7218008) q[2];
rz(-1.7828364) q[3];
sx q[3];
rz(-2.8163781) q[3];
sx q[3];
rz(1.1903919) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6841986) q[0];
sx q[0];
rz(-1.1379108) q[0];
sx q[0];
rz(-2.6275291) q[0];
rz(0.38814107) q[1];
sx q[1];
rz(-2.5297647) q[1];
sx q[1];
rz(-1.2302037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83004763) q[0];
sx q[0];
rz(-1.8732949) q[0];
sx q[0];
rz(1.3001469) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1934727) q[2];
sx q[2];
rz(-0.88250752) q[2];
sx q[2];
rz(1.4982868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2588604) q[1];
sx q[1];
rz(-0.89619918) q[1];
sx q[1];
rz(-2.0432977) q[1];
x q[2];
rz(2.4711126) q[3];
sx q[3];
rz(-1.3904212) q[3];
sx q[3];
rz(-0.84202858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6396883) q[2];
sx q[2];
rz(-1.014726) q[2];
sx q[2];
rz(0.91628966) q[2];
rz(-0.4942975) q[3];
sx q[3];
rz(-1.9435792) q[3];
sx q[3];
rz(-0.77427197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2275527) q[0];
sx q[0];
rz(-2.1844449) q[0];
sx q[0];
rz(-0.48602948) q[0];
rz(0.60274974) q[1];
sx q[1];
rz(-1.175368) q[1];
sx q[1];
rz(1.1154307) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19768342) q[0];
sx q[0];
rz(-0.75431529) q[0];
sx q[0];
rz(-0.59342845) q[0];
x q[1];
rz(0.40256315) q[2];
sx q[2];
rz(-2.3057115) q[2];
sx q[2];
rz(-2.1083567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.19957) q[1];
sx q[1];
rz(-1.1567819) q[1];
sx q[1];
rz(2.4041818) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0511469) q[3];
sx q[3];
rz(-1.8060038) q[3];
sx q[3];
rz(3.1344828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1439765) q[2];
sx q[2];
rz(-0.3207427) q[2];
sx q[2];
rz(2.0799267) q[2];
rz(1.5205421) q[3];
sx q[3];
rz(-1.1989667) q[3];
sx q[3];
rz(2.2475713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0735556) q[0];
sx q[0];
rz(-3.1017922) q[0];
sx q[0];
rz(1.918248) q[0];
rz(2.6262737) q[1];
sx q[1];
rz(-0.66910187) q[1];
sx q[1];
rz(1.3656778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66636234) q[0];
sx q[0];
rz(-1.7710553) q[0];
sx q[0];
rz(0.11563511) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3884945) q[2];
sx q[2];
rz(-1.4788091) q[2];
sx q[2];
rz(-2.0192041) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15348831) q[1];
sx q[1];
rz(-2.0447013) q[1];
sx q[1];
rz(-1.7405645) q[1];
rz(-pi) q[2];
rz(1.5346157) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.0042013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0627275) q[2];
sx q[2];
rz(-1.3389503) q[2];
sx q[2];
rz(1.8143181) q[2];
rz(1.7257388) q[3];
sx q[3];
rz(-0.073181987) q[3];
sx q[3];
rz(0.88254005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0075204) q[0];
sx q[0];
rz(-0.53711689) q[0];
sx q[0];
rz(2.9448331) q[0];
rz(0.21601954) q[1];
sx q[1];
rz(-1.0088423) q[1];
sx q[1];
rz(-2.3541727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8674054) q[0];
sx q[0];
rz(-1.3816773) q[0];
sx q[0];
rz(-1.8195137) q[0];
rz(2.6419053) q[2];
sx q[2];
rz(-2.5737107) q[2];
sx q[2];
rz(-0.5918006) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.710121) q[1];
sx q[1];
rz(-0.97968757) q[1];
sx q[1];
rz(1.5159056) q[1];
rz(-pi) q[2];
rz(1.2489572) q[3];
sx q[3];
rz(-1.3909381) q[3];
sx q[3];
rz(-1.091979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54062033) q[2];
sx q[2];
rz(-1.629402) q[2];
sx q[2];
rz(2.7579894) q[2];
rz(0.90841928) q[3];
sx q[3];
rz(-2.3288265) q[3];
sx q[3];
rz(-2.9878591) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28021321) q[0];
sx q[0];
rz(-0.52556831) q[0];
sx q[0];
rz(0.63661611) q[0];
rz(0.55066291) q[1];
sx q[1];
rz(-1.6267585) q[1];
sx q[1];
rz(-2.6689463) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.11926) q[0];
sx q[0];
rz(-0.68281931) q[0];
sx q[0];
rz(-0.14108087) q[0];
rz(-pi) q[1];
rz(3.0483629) q[2];
sx q[2];
rz(-1.5179253) q[2];
sx q[2];
rz(-0.44651647) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0090985) q[1];
sx q[1];
rz(-0.14202296) q[1];
sx q[1];
rz(0.36722398) q[1];
rz(2.5100559) q[3];
sx q[3];
rz(-1.8145921) q[3];
sx q[3];
rz(0.35758986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7134646) q[2];
sx q[2];
rz(-1.8992004) q[2];
sx q[2];
rz(3.0430651) q[2];
rz(0.030700961) q[3];
sx q[3];
rz(-1.5701598) q[3];
sx q[3];
rz(2.9876685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5310265) q[0];
sx q[0];
rz(-2.1838146) q[0];
sx q[0];
rz(2.4089693) q[0];
rz(3.1351807) q[1];
sx q[1];
rz(-2.1156204) q[1];
sx q[1];
rz(-1.7195255) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32860562) q[0];
sx q[0];
rz(-1.3805132) q[0];
sx q[0];
rz(1.3178131) q[0];
rz(-1.9344011) q[2];
sx q[2];
rz(-1.0997314) q[2];
sx q[2];
rz(0.94719145) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8274424) q[1];
sx q[1];
rz(-1.536751) q[1];
sx q[1];
rz(-2.9877547) q[1];
rz(-1.645817) q[3];
sx q[3];
rz(-1.0524233) q[3];
sx q[3];
rz(-3.1278267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0473359) q[2];
sx q[2];
rz(-1.7862659) q[2];
sx q[2];
rz(0.52200851) q[2];
rz(0.020307288) q[3];
sx q[3];
rz(-2.4419407) q[3];
sx q[3];
rz(-2.844753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.2688399) q[0];
sx q[0];
rz(-1.6167384) q[0];
sx q[0];
rz(-2.010345) q[0];
rz(2.3616683) q[1];
sx q[1];
rz(-1.2029388) q[1];
sx q[1];
rz(-2.6663229) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72544725) q[0];
sx q[0];
rz(-2.3088881) q[0];
sx q[0];
rz(1.0698143) q[0];
rz(-pi) q[1];
rz(-2.9067801) q[2];
sx q[2];
rz(-2.2249892) q[2];
sx q[2];
rz(-1.393911) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1058029) q[1];
sx q[1];
rz(-0.33638182) q[1];
sx q[1];
rz(2.2835284) q[1];
x q[2];
rz(1.4251208) q[3];
sx q[3];
rz(-2.2620107) q[3];
sx q[3];
rz(-2.447405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37107006) q[2];
sx q[2];
rz(-1.5757685) q[2];
sx q[2];
rz(-1.2062262) q[2];
rz(1.3953588) q[3];
sx q[3];
rz(-0.54308707) q[3];
sx q[3];
rz(3.0622845) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410011) q[0];
sx q[0];
rz(-0.19484367) q[0];
sx q[0];
rz(2.3051443) q[0];
rz(2.1438228) q[1];
sx q[1];
rz(-2.0496968) q[1];
sx q[1];
rz(-3.0539378) q[1];
rz(-0.069167698) q[2];
sx q[2];
rz(-1.3528878) q[2];
sx q[2];
rz(2.4895346) q[2];
rz(0.67881696) q[3];
sx q[3];
rz(-1.4842508) q[3];
sx q[3];
rz(2.1723464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
