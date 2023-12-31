OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(0.58340573) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(2.2489927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8164506) q[0];
sx q[0];
rz(-0.6439552) q[0];
sx q[0];
rz(-2.08026) q[0];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8436369) q[1];
sx q[1];
rz(-1.8814109) q[1];
sx q[1];
rz(-2.288504) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0953232) q[3];
sx q[3];
rz(-2.390063) q[3];
sx q[3];
rz(-0.7788333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2261752) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083667) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(0.57587409) q[0];
rz(-1.8946164) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.1670246) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562247) q[0];
sx q[0];
rz(-1.313374) q[0];
sx q[0];
rz(-0.060398922) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1738271) q[2];
sx q[2];
rz(-2.5644828) q[2];
sx q[2];
rz(1.4156262) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6845219) q[1];
sx q[1];
rz(-1.5091389) q[1];
sx q[1];
rz(0.80656273) q[1];
rz(2.5490709) q[3];
sx q[3];
rz(-0.7624818) q[3];
sx q[3];
rz(0.034686397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-2.8607821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4784933) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(-0.40027174) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(1.0148369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0815711) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(1.5006256) q[0];
x q[1];
rz(1.9760518) q[2];
sx q[2];
rz(-1.3477256) q[2];
sx q[2];
rz(0.16135339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9008357) q[1];
sx q[1];
rz(-1.6512617) q[1];
sx q[1];
rz(-2.1686173) q[1];
x q[2];
rz(-0.88481836) q[3];
sx q[3];
rz(-1.775145) q[3];
sx q[3];
rz(2.200978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(0.8992368) q[2];
rz(-2.4441161) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4863481) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-0.43193257) q[0];
rz(-0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-2.5057709) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0974554) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(1.9299279) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9767838) q[2];
sx q[2];
rz(-0.85068446) q[2];
sx q[2];
rz(1.009843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7446049) q[1];
sx q[1];
rz(-2.3481391) q[1];
sx q[1];
rz(0.26054392) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7838578) q[3];
sx q[3];
rz(-1.2902998) q[3];
sx q[3];
rz(-0.76663843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2137961) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(-0.33411807) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.14389811) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(-1.0744263) q[0];
rz(-0.74514666) q[1];
sx q[1];
rz(-1.4829758) q[1];
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
rz(-1.9946949) q[0];
sx q[0];
rz(-2.3197078) q[0];
rz(2.9859221) q[2];
sx q[2];
rz(-2.5241733) q[2];
sx q[2];
rz(1.9215259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3814195) q[1];
sx q[1];
rz(-0.96254327) q[1];
sx q[1];
rz(2.6890523) q[1];
x q[2];
rz(0.81551084) q[3];
sx q[3];
rz(-1.9933356) q[3];
sx q[3];
rz(-2.5583207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(0.29423514) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(0.0090573514) q[0];
rz(2.5065705) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(0.10805282) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51443433) q[0];
sx q[0];
rz(-1.2540199) q[0];
sx q[0];
rz(2.6677674) q[0];
x q[1];
rz(0.38185264) q[2];
sx q[2];
rz(-1.6994886) q[2];
sx q[2];
rz(-0.51634386) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5080155) q[1];
sx q[1];
rz(-2.0532236) q[1];
sx q[1];
rz(-1.1303933) q[1];
x q[2];
rz(1.990854) q[3];
sx q[3];
rz(-0.58745158) q[3];
sx q[3];
rz(2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(-0.078401119) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(-1.3972067) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-0.89362842) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2715627) q[0];
sx q[0];
rz(-1.8989519) q[0];
sx q[0];
rz(2.4782655) q[0];
rz(-pi) q[1];
rz(-1.0929843) q[2];
sx q[2];
rz(-1.3078948) q[2];
sx q[2];
rz(-0.36793567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9295846) q[1];
sx q[1];
rz(-1.7716904) q[1];
sx q[1];
rz(2.4561873) q[1];
rz(3.0393533) q[3];
sx q[3];
rz(-0.7623626) q[3];
sx q[3];
rz(1.8765212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(2.612109) q[2];
rz(-2.6654065) q[3];
sx q[3];
rz(-1.3323077) q[3];
sx q[3];
rz(0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(-1.09028) q[0];
rz(-0.11225637) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.9876678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.740828) q[0];
sx q[0];
rz(-0.65069288) q[0];
sx q[0];
rz(3.0853737) q[0];
rz(-pi) q[1];
rz(-1.7190785) q[2];
sx q[2];
rz(-1.0245171) q[2];
sx q[2];
rz(0.15020457) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0358419) q[1];
sx q[1];
rz(-0.81139794) q[1];
sx q[1];
rz(1.6559385) q[1];
rz(-pi) q[2];
rz(-2.047394) q[3];
sx q[3];
rz(-0.72924858) q[3];
sx q[3];
rz(-0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.551679) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(-0.38044688) q[2];
rz(2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8001051) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.170084) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78478783) q[0];
sx q[0];
rz(-2.6492282) q[0];
sx q[0];
rz(-0.82226336) q[0];
x q[1];
rz(1.3001928) q[2];
sx q[2];
rz(-2.3286208) q[2];
sx q[2];
rz(-1.9780618) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3635892) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(-2.7705454) q[1];
x q[2];
rz(1.3977259) q[3];
sx q[3];
rz(-0.67865463) q[3];
sx q[3];
rz(0.21733397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3141979) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(0.38273746) q[2];
rz(0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(0.25892192) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(-2.6616667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057248) q[0];
sx q[0];
rz(-0.94434443) q[0];
sx q[0];
rz(-0.41525526) q[0];
x q[1];
rz(-1.2725699) q[2];
sx q[2];
rz(-0.68694653) q[2];
sx q[2];
rz(-0.62703122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.568798) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(0.6154284) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2445573) q[3];
sx q[3];
rz(-1.9978791) q[3];
sx q[3];
rz(-2.2626812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2991128) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(-1.3170362) q[2];
rz(1.8995829) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(-2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
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
rz(2.6600044) q[2];
sx q[2];
rz(-1.7272186) q[2];
sx q[2];
rz(-2.061736) q[2];
rz(-1.8388207) q[3];
sx q[3];
rz(-1.839073) q[3];
sx q[3];
rz(3.0243235) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
