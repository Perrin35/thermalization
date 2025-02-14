OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2630513) q[0];
sx q[0];
rz(3.75293) q[0];
sx q[0];
rz(11.073025) q[0];
rz(0.85236323) q[1];
sx q[1];
rz(-1.3937997) q[1];
sx q[1];
rz(-1.2300904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2821101) q[0];
sx q[0];
rz(-1.9766269) q[0];
sx q[0];
rz(-0.52284436) q[0];
rz(-pi) q[1];
rz(1.1793433) q[2];
sx q[2];
rz(-1.9916996) q[2];
sx q[2];
rz(-0.98513033) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7802496) q[1];
sx q[1];
rz(-1.2318101) q[1];
sx q[1];
rz(-1.0913244) q[1];
x q[2];
rz(-1.8147179) q[3];
sx q[3];
rz(-2.8496242) q[3];
sx q[3];
rz(-1.0002182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7386231) q[2];
sx q[2];
rz(-1.3961184) q[2];
sx q[2];
rz(1.1633066) q[2];
rz(-0.93825424) q[3];
sx q[3];
rz(-2.3353751) q[3];
sx q[3];
rz(-1.7136542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52629483) q[0];
sx q[0];
rz(-1.2328923) q[0];
sx q[0];
rz(1.9799318) q[0];
rz(1.42234) q[1];
sx q[1];
rz(-1.6484345) q[1];
sx q[1];
rz(3.0551522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42936642) q[0];
sx q[0];
rz(-1.4131568) q[0];
sx q[0];
rz(-2.1033573) q[0];
rz(-pi) q[1];
rz(2.739846) q[2];
sx q[2];
rz(-0.52626393) q[2];
sx q[2];
rz(-1.9578735) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5878488) q[1];
sx q[1];
rz(-1.598473) q[1];
sx q[1];
rz(-1.5553655) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7419613) q[3];
sx q[3];
rz(-1.5054718) q[3];
sx q[3];
rz(0.32495503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5689759) q[2];
sx q[2];
rz(-0.47875753) q[2];
sx q[2];
rz(-3.0968481) q[2];
rz(2.4217126) q[3];
sx q[3];
rz(-1.6241112) q[3];
sx q[3];
rz(-1.450479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2429263) q[0];
sx q[0];
rz(-1.1761605) q[0];
sx q[0];
rz(1.2258919) q[0];
rz(-2.2236845) q[1];
sx q[1];
rz(-1.8501015) q[1];
sx q[1];
rz(-1.8570522) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51581406) q[0];
sx q[0];
rz(-0.5042133) q[0];
sx q[0];
rz(1.4674241) q[0];
x q[1];
rz(2.5199408) q[2];
sx q[2];
rz(-1.5772515) q[2];
sx q[2];
rz(2.9916951) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.94938147) q[1];
sx q[1];
rz(-0.95451372) q[1];
sx q[1];
rz(0.90009113) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35090943) q[3];
sx q[3];
rz(-2.4882462) q[3];
sx q[3];
rz(0.86337435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7104177) q[2];
sx q[2];
rz(-1.9937036) q[2];
sx q[2];
rz(2.5244024) q[2];
rz(2.2271473) q[3];
sx q[3];
rz(-2.4121598) q[3];
sx q[3];
rz(2.7845553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0451839) q[0];
sx q[0];
rz(-0.04627385) q[0];
sx q[0];
rz(2.1531877) q[0];
rz(1.660396) q[1];
sx q[1];
rz(-2.587187) q[1];
sx q[1];
rz(2.5344417) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0074127277) q[0];
sx q[0];
rz(-1.3152342) q[0];
sx q[0];
rz(-1.289053) q[0];
rz(-2.5012403) q[2];
sx q[2];
rz(-1.8812064) q[2];
sx q[2];
rz(-2.7057757) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1478575) q[1];
sx q[1];
rz(-2.5566935) q[1];
sx q[1];
rz(2.2123076) q[1];
x q[2];
rz(-1.0801717) q[3];
sx q[3];
rz(-2.9376224) q[3];
sx q[3];
rz(2.5734165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0561169) q[2];
sx q[2];
rz(-1.3885219) q[2];
sx q[2];
rz(0.2307387) q[2];
rz(0.76656109) q[3];
sx q[3];
rz(-0.14403382) q[3];
sx q[3];
rz(-1.2994331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6799927) q[0];
sx q[0];
rz(-1.5138641) q[0];
sx q[0];
rz(-3.0738714) q[0];
rz(1.8374775) q[1];
sx q[1];
rz(-1.8505406) q[1];
sx q[1];
rz(-1.7052604) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7505813) q[0];
sx q[0];
rz(-0.89871797) q[0];
sx q[0];
rz(2.6425022) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7737232) q[2];
sx q[2];
rz(-0.5022011) q[2];
sx q[2];
rz(-1.3853488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41714121) q[1];
sx q[1];
rz(-1.0994689) q[1];
sx q[1];
rz(-0.36496867) q[1];
rz(-pi) q[2];
rz(1.93238) q[3];
sx q[3];
rz(-1.7783594) q[3];
sx q[3];
rz(0.62011564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0962254) q[2];
sx q[2];
rz(-2.2787091) q[2];
sx q[2];
rz(-0.53675845) q[2];
rz(-1.2196352) q[3];
sx q[3];
rz(-0.9044956) q[3];
sx q[3];
rz(-1.9513244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9959975) q[0];
sx q[0];
rz(-2.1639316) q[0];
sx q[0];
rz(-1.4033432) q[0];
rz(-2.6373236) q[1];
sx q[1];
rz(-1.206617) q[1];
sx q[1];
rz(2.9020342) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1439622) q[0];
sx q[0];
rz(-0.24050783) q[0];
sx q[0];
rz(2.2872187) q[0];
rz(-pi) q[1];
rz(-2.1920565) q[2];
sx q[2];
rz(-2.4944948) q[2];
sx q[2];
rz(0.71352772) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18232432) q[1];
sx q[1];
rz(-2.2367918) q[1];
sx q[1];
rz(1.6823425) q[1];
rz(1.0658748) q[3];
sx q[3];
rz(-1.2137102) q[3];
sx q[3];
rz(0.034402196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2704894) q[2];
sx q[2];
rz(-1.5364545) q[2];
sx q[2];
rz(-2.8883873) q[2];
rz(-0.95868715) q[3];
sx q[3];
rz(-2.2264693) q[3];
sx q[3];
rz(-1.8668713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7764928) q[0];
sx q[0];
rz(-2.4204142) q[0];
sx q[0];
rz(-1.7505296) q[0];
rz(2.5513388) q[1];
sx q[1];
rz(-1.2140467) q[1];
sx q[1];
rz(-0.50277695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816526) q[0];
sx q[0];
rz(-1.556373) q[0];
sx q[0];
rz(-1.4653066) q[0];
rz(0.46765695) q[2];
sx q[2];
rz(-2.4291647) q[2];
sx q[2];
rz(-1.4234655) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5107837) q[1];
sx q[1];
rz(-1.8165339) q[1];
sx q[1];
rz(2.3151957) q[1];
rz(-pi) q[2];
x q[2];
rz(2.598113) q[3];
sx q[3];
rz(-1.8414652) q[3];
sx q[3];
rz(1.9668996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34660029) q[2];
sx q[2];
rz(-1.9575926) q[2];
sx q[2];
rz(-2.6896175) q[2];
rz(2.82708) q[3];
sx q[3];
rz(-1.9032685) q[3];
sx q[3];
rz(3.0934635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6314342) q[0];
sx q[0];
rz(-0.87961125) q[0];
sx q[0];
rz(-1.6360224) q[0];
rz(0.10558852) q[1];
sx q[1];
rz(-1.0786723) q[1];
sx q[1];
rz(2.4619335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15679793) q[0];
sx q[0];
rz(-2.121006) q[0];
sx q[0];
rz(-2.2722428) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7397268) q[2];
sx q[2];
rz(-0.99152126) q[2];
sx q[2];
rz(-1.4039672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9951156) q[1];
sx q[1];
rz(-1.5739643) q[1];
sx q[1];
rz(2.0407487) q[1];
rz(-pi) q[2];
rz(-1.7962436) q[3];
sx q[3];
rz(-1.7525008) q[3];
sx q[3];
rz(-0.34442201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52123657) q[2];
sx q[2];
rz(-2.607589) q[2];
sx q[2];
rz(2.3109069) q[2];
rz(1.0384809) q[3];
sx q[3];
rz(-2.4887648) q[3];
sx q[3];
rz(1.4000819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1164383) q[0];
sx q[0];
rz(-1.1517628) q[0];
sx q[0];
rz(1.0333767) q[0];
rz(2.0694464) q[1];
sx q[1];
rz(-1.1786345) q[1];
sx q[1];
rz(0.56195608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4406703) q[0];
sx q[0];
rz(-1.431968) q[0];
sx q[0];
rz(0.85900659) q[0];
x q[1];
rz(-0.91479723) q[2];
sx q[2];
rz(-1.2848228) q[2];
sx q[2];
rz(0.15483072) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1984673) q[1];
sx q[1];
rz(-2.2930305) q[1];
sx q[1];
rz(-2.6814382) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0086781541) q[3];
sx q[3];
rz(-1.5132705) q[3];
sx q[3];
rz(0.044139095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.098103913) q[2];
sx q[2];
rz(-1.4926814) q[2];
sx q[2];
rz(1.4081504) q[2];
rz(3.0812541) q[3];
sx q[3];
rz(-1.8294168) q[3];
sx q[3];
rz(-1.8797125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9982346) q[0];
sx q[0];
rz(-2.3682605) q[0];
sx q[0];
rz(1.3704569) q[0];
rz(-1.0073608) q[1];
sx q[1];
rz(-1.3887364) q[1];
sx q[1];
rz(-0.14152424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816982) q[0];
sx q[0];
rz(-2.0260677) q[0];
sx q[0];
rz(1.0819525) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3958489) q[2];
sx q[2];
rz(-0.088938449) q[2];
sx q[2];
rz(2.3740785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2685769) q[1];
sx q[1];
rz(-1.1339703) q[1];
sx q[1];
rz(2.718513) q[1];
rz(1.0753638) q[3];
sx q[3];
rz(-1.5902219) q[3];
sx q[3];
rz(2.7904449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1141899) q[2];
sx q[2];
rz(-1.9582615) q[2];
sx q[2];
rz(-1.2664504) q[2];
rz(-0.13628422) q[3];
sx q[3];
rz(-0.91629052) q[3];
sx q[3];
rz(-0.19792476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4895353) q[0];
sx q[0];
rz(-2.1283524) q[0];
sx q[0];
rz(1.5555489) q[0];
rz(2.1736705) q[1];
sx q[1];
rz(-0.5236917) q[1];
sx q[1];
rz(2.0655469) q[1];
rz(-1.684741) q[2];
sx q[2];
rz(-2.3782627) q[2];
sx q[2];
rz(-1.0397604) q[2];
rz(-1.8347673) q[3];
sx q[3];
rz(-0.19656678) q[3];
sx q[3];
rz(-1.8241775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
