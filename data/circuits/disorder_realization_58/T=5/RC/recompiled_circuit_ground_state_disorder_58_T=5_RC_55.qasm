OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4492884) q[0];
sx q[0];
rz(-2.6753354) q[0];
sx q[0];
rz(1.2944846) q[0];
rz(-0.26588384) q[1];
sx q[1];
rz(-2.9335913) q[1];
sx q[1];
rz(1.2535569) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3304015) q[0];
sx q[0];
rz(-1.6881144) q[0];
sx q[0];
rz(-0.53240029) q[0];
rz(-pi) q[1];
rz(0.8280379) q[2];
sx q[2];
rz(-2.5749505) q[2];
sx q[2];
rz(0.41929454) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2840062) q[1];
sx q[1];
rz(-2.2145971) q[1];
sx q[1];
rz(-2.8394152) q[1];
x q[2];
rz(-1.8951224) q[3];
sx q[3];
rz(-1.5889349) q[3];
sx q[3];
rz(-0.42551009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1050538) q[2];
sx q[2];
rz(-1.7221071) q[2];
sx q[2];
rz(1.4744021) q[2];
rz(0.27753943) q[3];
sx q[3];
rz(-1.8614635) q[3];
sx q[3];
rz(2.2138219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14811806) q[0];
sx q[0];
rz(-0.19522218) q[0];
sx q[0];
rz(1.8147722) q[0];
rz(-2.7534292) q[1];
sx q[1];
rz(-1.6366942) q[1];
sx q[1];
rz(0.51663748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46088567) q[0];
sx q[0];
rz(-0.26038489) q[0];
sx q[0];
rz(2.5234114) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1519127) q[2];
sx q[2];
rz(-1.7563213) q[2];
sx q[2];
rz(-2.1806661) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7611446) q[1];
sx q[1];
rz(-0.7711322) q[1];
sx q[1];
rz(1.3517595) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2907007) q[3];
sx q[3];
rz(-1.3317654) q[3];
sx q[3];
rz(1.7861507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8119729) q[2];
sx q[2];
rz(-2.6361578) q[2];
sx q[2];
rz(2.2774515) q[2];
rz(2.4498074) q[3];
sx q[3];
rz(-2.0103879) q[3];
sx q[3];
rz(-2.2085021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8253887) q[0];
sx q[0];
rz(-0.80556691) q[0];
sx q[0];
rz(1.7387996) q[0];
rz(0.56651506) q[1];
sx q[1];
rz(-1.6367876) q[1];
sx q[1];
rz(-1.8623955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7544781) q[0];
sx q[0];
rz(-0.25969782) q[0];
sx q[0];
rz(-1.5952871) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6635936) q[2];
sx q[2];
rz(-1.0339811) q[2];
sx q[2];
rz(-1.9271242) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6237026) q[1];
sx q[1];
rz(-1.5346171) q[1];
sx q[1];
rz(2.0149798) q[1];
rz(-pi) q[2];
rz(-2.6353929) q[3];
sx q[3];
rz(-2.4512614) q[3];
sx q[3];
rz(1.797628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0731571) q[2];
sx q[2];
rz(-2.7918039) q[2];
sx q[2];
rz(-1.5901828) q[2];
rz(0.50897151) q[3];
sx q[3];
rz(-0.83566982) q[3];
sx q[3];
rz(1.6905748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1125672) q[0];
sx q[0];
rz(-1.4968766) q[0];
sx q[0];
rz(-0.56370869) q[0];
rz(1.6911223) q[1];
sx q[1];
rz(-1.1553973) q[1];
sx q[1];
rz(-2.3203826) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9552069) q[0];
sx q[0];
rz(-1.5599376) q[0];
sx q[0];
rz(-0.71594724) q[0];
x q[1];
rz(-2.8217373) q[2];
sx q[2];
rz(-1.8531688) q[2];
sx q[2];
rz(-0.47114633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5933696) q[1];
sx q[1];
rz(-2.0851622) q[1];
sx q[1];
rz(-2.3002808) q[1];
rz(-pi) q[2];
rz(-1.508237) q[3];
sx q[3];
rz(-0.54397115) q[3];
sx q[3];
rz(-2.6450752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.700909) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(1.7636501) q[2];
rz(2.9772229) q[3];
sx q[3];
rz(-1.5956968) q[3];
sx q[3];
rz(-0.36851287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57581562) q[0];
sx q[0];
rz(-1.4482647) q[0];
sx q[0];
rz(2.6337295) q[0];
rz(0.85743088) q[1];
sx q[1];
rz(-1.6897759) q[1];
sx q[1];
rz(-0.47971496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21151152) q[0];
sx q[0];
rz(-0.9418315) q[0];
sx q[0];
rz(-1.1674561) q[0];
rz(-pi) q[1];
rz(1.4998593) q[2];
sx q[2];
rz(-1.1878769) q[2];
sx q[2];
rz(1.4986582) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.060412571) q[1];
sx q[1];
rz(-0.65250373) q[1];
sx q[1];
rz(-1.9461701) q[1];
rz(-pi) q[2];
rz(-0.028939441) q[3];
sx q[3];
rz(-1.2914198) q[3];
sx q[3];
rz(-2.8795163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63626426) q[2];
sx q[2];
rz(-2.0365066) q[2];
sx q[2];
rz(0.21978933) q[2];
rz(-2.9124741) q[3];
sx q[3];
rz(-0.90355211) q[3];
sx q[3];
rz(-0.98178274) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.583928) q[0];
sx q[0];
rz(-2.5983577) q[0];
sx q[0];
rz(1.1274717) q[0];
rz(-2.8358031) q[1];
sx q[1];
rz(-1.1916279) q[1];
sx q[1];
rz(2.9514899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7641417) q[0];
sx q[0];
rz(-0.5849291) q[0];
sx q[0];
rz(-2.1795953) q[0];
rz(-pi) q[1];
rz(-2.7626412) q[2];
sx q[2];
rz(-2.3109461) q[2];
sx q[2];
rz(-1.3791858) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1446962) q[1];
sx q[1];
rz(-2.3544925) q[1];
sx q[1];
rz(1.1622278) q[1];
x q[2];
rz(-0.28921952) q[3];
sx q[3];
rz(-0.26258024) q[3];
sx q[3];
rz(-2.5724831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1216792) q[2];
sx q[2];
rz(-1.6807669) q[2];
sx q[2];
rz(-1.100568) q[2];
rz(1.8371643) q[3];
sx q[3];
rz(-1.5521939) q[3];
sx q[3];
rz(-1.7884375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28534999) q[0];
sx q[0];
rz(-2.7142363) q[0];
sx q[0];
rz(2.5653978) q[0];
rz(1.0528437) q[1];
sx q[1];
rz(-2.5710227) q[1];
sx q[1];
rz(-1.0594692) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16082009) q[0];
sx q[0];
rz(-2.4376025) q[0];
sx q[0];
rz(-2.9454548) q[0];
x q[1];
rz(3.1131831) q[2];
sx q[2];
rz(-1.1724262) q[2];
sx q[2];
rz(2.0223126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1394315) q[1];
sx q[1];
rz(-2.0116848) q[1];
sx q[1];
rz(-1.7825761) q[1];
x q[2];
rz(0.93848159) q[3];
sx q[3];
rz(-0.96358591) q[3];
sx q[3];
rz(-0.069165088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2759555) q[2];
sx q[2];
rz(-2.553678) q[2];
sx q[2];
rz(2.0666583) q[2];
rz(0.86447191) q[3];
sx q[3];
rz(-1.8755707) q[3];
sx q[3];
rz(1.5786494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39701715) q[0];
sx q[0];
rz(-0.96857849) q[0];
sx q[0];
rz(-2.6343935) q[0];
rz(-1.1811258) q[1];
sx q[1];
rz(-1.6543417) q[1];
sx q[1];
rz(-2.2307253) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1169029) q[0];
sx q[0];
rz(-2.0729613) q[0];
sx q[0];
rz(-0.85310081) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8751926) q[2];
sx q[2];
rz(-2.5385529) q[2];
sx q[2];
rz(1.837544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6421318) q[1];
sx q[1];
rz(-1.5593464) q[1];
sx q[1];
rz(1.4169372) q[1];
x q[2];
rz(-2.8975119) q[3];
sx q[3];
rz(-1.968813) q[3];
sx q[3];
rz(0.73792968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38535038) q[2];
sx q[2];
rz(-1.664398) q[2];
sx q[2];
rz(0.64794668) q[2];
rz(3.044965) q[3];
sx q[3];
rz(-2.1164618) q[3];
sx q[3];
rz(2.1881762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18428093) q[0];
sx q[0];
rz(-2.1537557) q[0];
sx q[0];
rz(1.3405569) q[0];
rz(0.52349177) q[1];
sx q[1];
rz(-1.610787) q[1];
sx q[1];
rz(-1.2045822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25239326) q[0];
sx q[0];
rz(-2.3601818) q[0];
sx q[0];
rz(-1.8535421) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7707945) q[2];
sx q[2];
rz(-1.424779) q[2];
sx q[2];
rz(-2.4869362) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2661926) q[1];
sx q[1];
rz(-1.0169694) q[1];
sx q[1];
rz(-1.758105) q[1];
rz(-pi) q[2];
rz(-0.91851652) q[3];
sx q[3];
rz(-0.61448408) q[3];
sx q[3];
rz(2.8826734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.039310731) q[2];
sx q[2];
rz(-1.0264341) q[2];
sx q[2];
rz(-0.2956051) q[2];
rz(3.0292656) q[3];
sx q[3];
rz(-1.5528409) q[3];
sx q[3];
rz(2.2789392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2304147) q[0];
sx q[0];
rz(-2.6417612) q[0];
sx q[0];
rz(-2.3590132) q[0];
rz(0.29620194) q[1];
sx q[1];
rz(-1.9007416) q[1];
sx q[1];
rz(-2.9414419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5458465) q[0];
sx q[0];
rz(-2.3538618) q[0];
sx q[0];
rz(2.911522) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28130071) q[2];
sx q[2];
rz(-2.1185115) q[2];
sx q[2];
rz(2.3980464) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0240062) q[1];
sx q[1];
rz(-0.96029687) q[1];
sx q[1];
rz(2.4897051) q[1];
rz(-2.5298821) q[3];
sx q[3];
rz(-1.0308305) q[3];
sx q[3];
rz(1.5276791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.079411) q[2];
sx q[2];
rz(-0.14180413) q[2];
sx q[2];
rz(1.0395435) q[2];
rz(-3.1145575) q[3];
sx q[3];
rz(-0.98374933) q[3];
sx q[3];
rz(-2.1496617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(2.1886002) q[0];
sx q[0];
rz(-2.47692) q[0];
sx q[0];
rz(1.1823786) q[0];
rz(1.4324808) q[1];
sx q[1];
rz(-1.2624337) q[1];
sx q[1];
rz(1.591325) q[1];
rz(2.9903632) q[2];
sx q[2];
rz(-1.4561903) q[2];
sx q[2];
rz(-1.0753808) q[2];
rz(-1.6662206) q[3];
sx q[3];
rz(-2.0153042) q[3];
sx q[3];
rz(2.5074625) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
