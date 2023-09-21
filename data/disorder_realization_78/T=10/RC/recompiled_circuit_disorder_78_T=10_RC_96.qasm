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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8164506) q[0];
sx q[0];
rz(-2.4976375) q[0];
sx q[0];
rz(-2.08026) q[0];
rz(-pi) q[1];
rz(-0.61172723) q[2];
sx q[2];
rz(-2.3731542) q[2];
sx q[2];
rz(-0.4380463) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6751911) q[1];
sx q[1];
rz(-2.2474504) q[1];
sx q[1];
rz(-2.7387709) q[1];
rz(-0.046269429) q[3];
sx q[3];
rz(-2.390063) q[3];
sx q[3];
rz(2.3627594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.1606476) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.083667) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(-2.5657186) q[0];
rz(1.2469762) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(1.974568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052159667) q[0];
sx q[0];
rz(-0.26425996) q[0];
sx q[0];
rz(1.7961851) q[0];
rz(2.1115233) q[2];
sx q[2];
rz(-1.7833372) q[2];
sx q[2];
rz(-2.9589047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0547304) q[1];
sx q[1];
rz(-2.3332101) q[1];
sx q[1];
rz(-3.0562835) q[1];
rz(-2.0608276) q[3];
sx q[3];
rz(-0.96066517) q[3];
sx q[3];
rz(-0.78435635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.033826753) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(3.068148) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.6630994) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(-2.1267557) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9432362) q[0];
sx q[0];
rz(-2.495932) q[0];
sx q[0];
rz(-0.093430324) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24200183) q[2];
sx q[2];
rz(-1.9654462) q[2];
sx q[2];
rz(1.5040656) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9008357) q[1];
sx q[1];
rz(-1.490331) q[1];
sx q[1];
rz(-2.1686173) q[1];
rz(-pi) q[2];
rz(2.879911) q[3];
sx q[3];
rz(-2.2398584) q[3];
sx q[3];
rz(0.79479587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0959452) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(2.2423559) q[2];
rz(0.69747654) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(-2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4863481) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(2.7096601) q[0];
rz(0.63255429) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(-2.5057709) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0441372) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(-1.2116648) q[0];
rz(0.16480883) q[2];
sx q[2];
rz(-2.2909082) q[2];
sx q[2];
rz(2.1317496) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7446049) q[1];
sx q[1];
rz(-0.79345353) q[1];
sx q[1];
rz(2.8810487) q[1];
rz(-pi) q[2];
rz(-1.7838578) q[3];
sx q[3];
rz(-1.2902998) q[3];
sx q[3];
rz(-2.3749542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9277966) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(2.4528743) q[2];
rz(0.33411807) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(-0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9976945) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(1.0744263) q[0];
rz(-2.396446) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(2.863046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700096) q[0];
sx q[0];
rz(-0.9013114) q[0];
sx q[0];
rz(2.5894126) q[0];
rz(-pi) q[1];
rz(2.9859221) q[2];
sx q[2];
rz(-0.61741932) q[2];
sx q[2];
rz(1.2200668) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7601732) q[1];
sx q[1];
rz(-2.1790494) q[1];
sx q[1];
rz(0.45254032) q[1];
rz(0.99028011) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(2.5648404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(-0.29423514) q[2];
rz(-3.0596628) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(3.1325353) q[0];
rz(2.5065705) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(-3.0335398) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2436284) q[0];
sx q[0];
rz(-2.0192696) q[0];
sx q[0];
rz(-1.2178221) q[0];
rz(-pi) q[1];
rz(0.33424218) q[2];
sx q[2];
rz(-0.40194449) q[2];
sx q[2];
rz(2.3964756) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.84032816) q[1];
sx q[1];
rz(-2.5003308) q[1];
sx q[1];
rz(2.458359) q[1];
x q[2];
rz(-1.990854) q[3];
sx q[3];
rz(-2.5541411) q[3];
sx q[3];
rz(2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(-1.8704869) q[2];
rz(-0.078401119) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(-2.0097849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8909797) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(2.4839731) q[0];
rz(-1.744386) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(2.2479642) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0890869) q[0];
sx q[0];
rz(-2.1930709) q[0];
sx q[0];
rz(-1.9786579) q[0];
rz(-1.0412752) q[2];
sx q[2];
rz(-0.5404226) q[2];
sx q[2];
rz(1.668001) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.212008) q[1];
sx q[1];
rz(-1.3699023) q[1];
sx q[1];
rz(-2.4561873) q[1];
x q[2];
rz(0.75974792) q[3];
sx q[3];
rz(-1.6413416) q[3];
sx q[3];
rz(-2.7618046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39067337) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(-0.52948362) q[2];
rz(-2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(2.7990394) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.9876678) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40076462) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(3.0853737) q[0];
rz(2.9032193) q[2];
sx q[2];
rz(-0.56406883) q[2];
sx q[2];
rz(2.7114045) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1057508) q[1];
sx q[1];
rz(-2.3301947) q[1];
sx q[1];
rz(1.4856542) q[1];
rz(-0.89973255) q[3];
sx q[3];
rz(-1.2601488) q[3];
sx q[3];
rz(-1.8048546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5899137) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(2.5019116) q[0];
rz(-1.9027963) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(1.9715086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4719452) q[0];
sx q[0];
rz(-1.8983316) q[0];
sx q[0];
rz(-1.1963084) q[0];
x q[1];
rz(-0.27530382) q[2];
sx q[2];
rz(-0.79553662) q[2];
sx q[2];
rz(-2.3616633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77800345) q[1];
sx q[1];
rz(-2.6215141) q[1];
sx q[1];
rz(2.7705454) q[1];
rz(-1.3977259) q[3];
sx q[3];
rz(-2.462938) q[3];
sx q[3];
rz(0.21733397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8273948) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-2.8826707) q[0];
rz(2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(0.47992596) q[1];
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
rz(-0.9058814) q[2];
sx q[2];
rz(-1.3833628) q[2];
sx q[2];
rz(0.71042176) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.568798) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(-2.5261643) q[1];
rz(-pi) q[2];
rz(0.89703538) q[3];
sx q[3];
rz(-1.1437136) q[3];
sx q[3];
rz(0.87891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(-2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9823572) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(2.6600044) q[2];
sx q[2];
rz(-1.7272186) q[2];
sx q[2];
rz(-2.061736) q[2];
rz(-2.8638774) q[3];
sx q[3];
rz(-1.82901) q[3];
sx q[3];
rz(1.3808586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
