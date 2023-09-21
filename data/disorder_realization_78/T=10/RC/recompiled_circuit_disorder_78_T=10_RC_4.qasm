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
rz(-2.5506033) q[0];
sx q[0];
rz(-0.58340573) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(2.2489927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8569782) q[0];
sx q[0];
rz(-1.0191139) q[0];
sx q[0];
rz(2.7906228) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66944389) q[2];
sx q[2];
rz(-1.9813683) q[2];
sx q[2];
rz(-2.4759811) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0639227) q[1];
sx q[1];
rz(-0.7709255) q[1];
sx q[1];
rz(-2.0248807) q[1];
rz(-0.046269429) q[3];
sx q[3];
rz(-0.75152961) q[3];
sx q[3];
rz(-2.3627594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9154174) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(-2.9246269) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0579257) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(-0.57587409) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(-1.1670246) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052159667) q[0];
sx q[0];
rz(-2.8773327) q[0];
sx q[0];
rz(1.7961851) q[0];
rz(-2.8950047) q[2];
sx q[2];
rz(-1.0435259) q[2];
sx q[2];
rz(-1.2621244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.963672) q[1];
sx q[1];
rz(-2.3753787) q[1];
sx q[1];
rz(-1.6597762) q[1];
rz(-pi) q[2];
rz(2.4715273) q[3];
sx q[3];
rz(-1.1747922) q[3];
sx q[3];
rz(2.0585287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(0.21437422) q[2];
rz(-3.068148) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(2.8607821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4784933) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.2269155) q[0];
rz(2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(2.1267557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0815711) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(-1.5006256) q[0];
x q[1];
rz(2.8995908) q[2];
sx q[2];
rz(-1.1761464) q[2];
sx q[2];
rz(1.637527) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6940569) q[1];
sx q[1];
rz(-0.60255614) q[1];
sx q[1];
rz(-1.7130997) q[1];
x q[2];
rz(0.26168163) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(-2.3467968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0456475) q[2];
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
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863481) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(2.7096601) q[0];
rz(-2.5090384) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(2.5057709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0441372) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(1.2116648) q[0];
rz(2.2976774) q[2];
sx q[2];
rz(-1.6944431) q[2];
sx q[2];
rz(2.4713949) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1078474) q[1];
sx q[1];
rz(-0.81106942) q[1];
sx q[1];
rz(1.8268405) q[1];
x q[2];
rz(-2.5084247) q[3];
sx q[3];
rz(-2.79106) q[3];
sx q[3];
rz(-1.4299973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9277966) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9976945) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(2.396446) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(2.863046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644088) q[0];
sx q[0];
rz(-1.9946949) q[0];
sx q[0];
rz(0.82188481) q[0];
x q[1];
rz(-2.5299046) q[2];
sx q[2];
rz(-1.6606765) q[2];
sx q[2];
rz(-0.22345605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6812233) q[1];
sx q[1];
rz(-1.2037828) q[1];
sx q[1];
rz(0.91194921) q[1];
rz(-2.1513125) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(-0.57675225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4429861) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(-2.8473575) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(3.1325353) q[0];
rz(-0.63502216) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(-0.10805282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8979643) q[0];
sx q[0];
rz(-2.0192696) q[0];
sx q[0];
rz(1.2178221) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7093541) q[2];
sx q[2];
rz(-1.9493305) q[2];
sx q[2];
rz(2.0356503) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6335771) q[1];
sx q[1];
rz(-2.0532236) q[1];
sx q[1];
rz(2.0111994) q[1];
x q[2];
rz(1.1507387) q[3];
sx q[3];
rz(-2.5541411) q[3];
sx q[3];
rz(-0.95668876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(-1.2711058) q[2];
rz(-0.078401119) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(0.65761956) q[0];
rz(1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-0.89362842) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0890869) q[0];
sx q[0];
rz(-2.1930709) q[0];
sx q[0];
rz(-1.1629348) q[0];
rz(-pi) q[1];
rz(2.0486084) q[2];
sx q[2];
rz(-1.8336979) q[2];
sx q[2];
rz(-2.773657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11941714) q[1];
sx q[1];
rz(-0.70964538) q[1];
sx q[1];
rz(2.8303353) q[1];
rz(-pi) q[2];
rz(-2.3818447) q[3];
sx q[3];
rz(-1.6413416) q[3];
sx q[3];
rz(0.37978803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39067337) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(2.612109) q[2];
rz(0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(-0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.7628409) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(0.11225637) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.1539248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40076462) q[0];
sx q[0];
rz(-0.65069288) q[0];
sx q[0];
rz(-3.0853737) q[0];
rz(-1.4225142) q[2];
sx q[2];
rz(-2.1170756) q[2];
sx q[2];
rz(-2.9913881) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1057508) q[1];
sx q[1];
rz(-2.3301947) q[1];
sx q[1];
rz(1.6559385) q[1];
x q[2];
rz(0.89973255) q[3];
sx q[3];
rz(-1.8814439) q[3];
sx q[3];
rz(-1.8048546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.551679) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(-0.38044688) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34148759) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(-1.9027963) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.9715086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1146678) q[0];
sx q[0];
rz(-1.2171193) q[0];
sx q[0];
rz(-0.35004079) q[0];
rz(-0.77634546) q[2];
sx q[2];
rz(-1.3753969) q[2];
sx q[2];
rz(-0.59567829) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6744086) q[1];
sx q[1];
rz(-1.7519752) q[1];
sx q[1];
rz(2.6513537) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3977259) q[3];
sx q[3];
rz(-2.462938) q[3];
sx q[3];
rz(2.9242587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8273948) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(2.2132204) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(0.66463566) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2255573) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(2.8826707) q[0];
rz(2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(0.47992596) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5878764) q[0];
sx q[0];
rz(-1.237861) q[0];
sx q[0];
rz(-2.2399708) q[0];
rz(-pi) q[1];
rz(-0.9058814) q[2];
sx q[2];
rz(-1.7582298) q[2];
sx q[2];
rz(2.4311709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2952134) q[1];
sx q[1];
rz(-2.0800253) q[1];
sx q[1];
rz(1.9766115) q[1];
rz(-pi) q[2];
rz(-2.2009833) q[3];
sx q[3];
rz(-2.3621231) q[3];
sx q[3];
rz(1.1704695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(-1.8995829) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15923545) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(1.7469035) q[2];
sx q[2];
rz(-2.0460143) q[2];
sx q[2];
rz(2.569414) q[2];
rz(1.3027719) q[3];
sx q[3];
rz(-1.839073) q[3];
sx q[3];
rz(3.0243235) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
