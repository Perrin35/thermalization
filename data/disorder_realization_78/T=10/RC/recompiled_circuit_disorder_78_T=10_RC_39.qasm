OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(3.732582) q[0];
sx q[0];
rz(8.8413722) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(-0.89259994) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8164506) q[0];
sx q[0];
rz(-0.6439552) q[0];
sx q[0];
rz(-1.0613326) q[0];
rz(-pi) q[1];
rz(-2.4721488) q[2];
sx q[2];
rz(-1.1602243) q[2];
sx q[2];
rz(-2.4759811) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6751911) q[1];
sx q[1];
rz(-2.2474504) q[1];
sx q[1];
rz(-2.7387709) q[1];
rz(-pi) q[2];
rz(3.0953232) q[3];
sx q[3];
rz(-2.390063) q[3];
sx q[3];
rz(-0.7788333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083667) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(-2.5657186) q[0];
rz(-1.2469762) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.974568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28536797) q[0];
sx q[0];
rz(-1.8282187) q[0];
sx q[0];
rz(3.0811937) q[0];
rz(-pi) q[1];
rz(-0.24658792) q[2];
sx q[2];
rz(-2.0980667) q[2];
sx q[2];
rz(1.8794683) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45707073) q[1];
sx q[1];
rz(-1.5091389) q[1];
sx q[1];
rz(-2.3350299) q[1];
rz(1.0807651) q[3];
sx q[3];
rz(-0.96066517) q[3];
sx q[3];
rz(-0.78435635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.033826753) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-0.21437422) q[2];
rz(-0.073444627) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.9146772) q[0];
rz(0.40027174) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(-1.0148369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0600216) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(1.5006256) q[0];
x q[1];
rz(0.24200183) q[2];
sx q[2];
rz(-1.1761464) q[2];
sx q[2];
rz(-1.637527) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6940569) q[1];
sx q[1];
rz(-2.5390365) q[1];
sx q[1];
rz(-1.7130997) q[1];
rz(-2.2567743) q[3];
sx q[3];
rz(-1.3664477) q[3];
sx q[3];
rz(2.200978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0456475) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(2.2423559) q[2];
rz(-2.4441161) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(2.7096601) q[0];
rz(-0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-2.5057709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5084002) q[0];
sx q[0];
rz(-2.0661372) q[0];
sx q[0];
rz(0.20423996) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84391528) q[2];
sx q[2];
rz(-1.4471495) q[2];
sx q[2];
rz(-2.4713949) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0337452) q[1];
sx q[1];
rz(-0.81106942) q[1];
sx q[1];
rz(-1.8268405) q[1];
rz(-0.28663978) q[3];
sx q[3];
rz(-1.7754103) q[3];
sx q[3];
rz(2.3972547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9277966) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(0.68871838) q[2];
rz(-2.8074746) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(-0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14389811) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(2.0671663) q[0];
rz(0.74514666) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(0.27854663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4771839) q[0];
sx q[0];
rz(-1.1468977) q[0];
sx q[0];
rz(2.3197078) q[0];
rz(1.4611545) q[2];
sx q[2];
rz(-0.96193681) q[2];
sx q[2];
rz(-1.4102175) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3814195) q[1];
sx q[1];
rz(-0.96254327) q[1];
sx q[1];
rz(2.6890523) q[1];
rz(-pi) q[2];
rz(2.1513125) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(0.57675225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(-0.29423514) q[2];
rz(3.0596628) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(-3.1053655) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1086403) q[0];
sx q[0];
rz(-2.2886798) q[0];
sx q[0];
rz(-3.1325353) q[0];
rz(-0.63502216) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(3.0335398) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5393339) q[0];
sx q[0];
rz(-0.56319153) q[0];
sx q[0];
rz(-2.5186033) q[0];
x q[1];
rz(-2.8073505) q[2];
sx q[2];
rz(-0.40194449) q[2];
sx q[2];
rz(-0.7451171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15247004) q[1];
sx q[1];
rz(-1.9580541) q[1];
sx q[1];
rz(-0.5247922) q[1];
x q[2];
rz(2.8764358) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(1.6924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(-0.078401119) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(2.4839731) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(-0.89362842) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0890869) q[0];
sx q[0];
rz(-2.1930709) q[0];
sx q[0];
rz(-1.9786579) q[0];
rz(-pi) q[1];
rz(-2.8473179) q[2];
sx q[2];
rz(-2.0308959) q[2];
sx q[2];
rz(1.0690881) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11941714) q[1];
sx q[1];
rz(-2.4319473) q[1];
sx q[1];
rz(-0.31125734) q[1];
x q[2];
rz(-3.0393533) q[3];
sx q[3];
rz(-0.7623626) q[3];
sx q[3];
rz(1.2650714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7509193) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(-2.612109) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3787518) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(1.09028) q[0];
rz(-3.0293363) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.1539248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9268122) q[0];
sx q[0];
rz(-1.5367537) q[0];
sx q[0];
rz(-0.649931) q[0];
rz(-0.55118982) q[2];
sx q[2];
rz(-1.6973719) q[2];
sx q[2];
rz(1.7984496) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47626074) q[1];
sx q[1];
rz(-1.6325103) q[1];
sx q[1];
rz(-0.76121059) q[1];
rz(-pi) q[2];
x q[2];
rz(2.047394) q[3];
sx q[3];
rz(-2.4123441) q[3];
sx q[3];
rz(3.0081528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-0.38044688) q[2];
rz(2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(0.63968101) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(1.9715086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78478783) q[0];
sx q[0];
rz(-0.49236449) q[0];
sx q[0];
rz(-0.82226336) q[0];
rz(1.8413999) q[2];
sx q[2];
rz(-2.3286208) q[2];
sx q[2];
rz(-1.1635309) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3635892) q[1];
sx q[1];
rz(-2.6215141) q[1];
sx q[1];
rz(-0.37104721) q[1];
x q[2];
rz(-2.2421078) q[3];
sx q[3];
rz(-1.4624819) q[3];
sx q[3];
rz(1.9233821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3141979) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(-0.38273746) q[2];
rz(-2.2132204) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(2.8826707) q[0];
rz(2.4312773) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(-0.47992596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5878764) q[0];
sx q[0];
rz(-1.237861) q[0];
sx q[0];
rz(-0.90162189) q[0];
rz(-2.9051022) q[2];
sx q[2];
rz(-2.2220526) q[2];
sx q[2];
rz(2.1361534) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2952134) q[1];
sx q[1];
rz(-2.0800253) q[1];
sx q[1];
rz(1.9766115) q[1];
rz(-pi) q[2];
rz(-2.6142526) q[3];
sx q[3];
rz(-2.1747327) q[3];
sx q[3];
rz(0.37249836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84247983) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(-1.8245565) q[2];
rz(1.2420098) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(-2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823572) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(-0.48158823) q[2];
sx q[2];
rz(-1.7272186) q[2];
sx q[2];
rz(-2.061736) q[2];
rz(-1.3027719) q[3];
sx q[3];
rz(-1.3025197) q[3];
sx q[3];
rz(-0.11726911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];