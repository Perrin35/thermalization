OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(6.364967) q[0];
sx q[0];
rz(9.9262417) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0599521) q[0];
sx q[0];
rz(-1.2498706) q[0];
sx q[0];
rz(3.0517464) q[0];
x q[1];
rz(-0.85876043) q[2];
sx q[2];
rz(-2.3220064) q[2];
sx q[2];
rz(-0.31529266) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8391708) q[1];
sx q[1];
rz(-1.2245373) q[1];
sx q[1];
rz(0.94567169) q[1];
x q[2];
rz(-0.97500719) q[3];
sx q[3];
rz(-2.4132204) q[3];
sx q[3];
rz(-1.7155852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(-0.55603975) q[2];
rz(-2.3089144) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(0.15727501) q[0];
rz(-2.8804624) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(0.10903407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5058274) q[0];
sx q[0];
rz(-0.286239) q[0];
sx q[0];
rz(2.2155227) q[0];
x q[1];
rz(-2.878506) q[2];
sx q[2];
rz(-1.3795985) q[2];
sx q[2];
rz(2.4441602) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82169689) q[1];
sx q[1];
rz(-0.79274717) q[1];
sx q[1];
rz(2.4577623) q[1];
rz(-pi) q[2];
rz(1.3974959) q[3];
sx q[3];
rz(-2.0850075) q[3];
sx q[3];
rz(-1.3483931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(1.8781352) q[2];
rz(-0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(-2.893977) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-0.48167357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047086296) q[0];
sx q[0];
rz(-1.9425834) q[0];
sx q[0];
rz(-1.0466172) q[0];
x q[1];
rz(-1.1142251) q[2];
sx q[2];
rz(-2.4764428) q[2];
sx q[2];
rz(-2.3014724) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86320089) q[1];
sx q[1];
rz(-0.87676261) q[1];
sx q[1];
rz(0.73514003) q[1];
rz(-pi) q[2];
rz(0.84677245) q[3];
sx q[3];
rz(-2.0623042) q[3];
sx q[3];
rz(1.5566952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3383011) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(0.49989191) q[2];
rz(0.56097427) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19701476) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(-0.552185) q[0];
rz(1.5532956) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.2447371) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9024076) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(-0.0048588077) q[0];
x q[1];
rz(-0.130529) q[2];
sx q[2];
rz(-1.7703238) q[2];
sx q[2];
rz(-0.56632698) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0075477608) q[1];
sx q[1];
rz(-1.782685) q[1];
sx q[1];
rz(-1.0210387) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7961058) q[3];
sx q[3];
rz(-2.1300737) q[3];
sx q[3];
rz(2.6252281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(-1.4771279) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(-0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(-3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(-1.5030456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3969288) q[0];
sx q[0];
rz(-1.443112) q[0];
sx q[0];
rz(-0.98608195) q[0];
x q[1];
rz(-2.224515) q[2];
sx q[2];
rz(-3.0055025) q[2];
sx q[2];
rz(1.9249141) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0283708) q[1];
sx q[1];
rz(-1.2591259) q[1];
sx q[1];
rz(0.10722864) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6081229) q[3];
sx q[3];
rz(-0.57538486) q[3];
sx q[3];
rz(1.433978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6333255) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(-1.903669) q[2];
rz(-1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(-0.52156633) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313107) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(2.8748728) q[0];
rz(2.5807014) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(0.7985324) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.610299) q[0];
sx q[0];
rz(-1.4302505) q[0];
sx q[0];
rz(1.2929686) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59136765) q[2];
sx q[2];
rz(-2.5632576) q[2];
sx q[2];
rz(-2.8665198) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9763899) q[1];
sx q[1];
rz(-2.5678647) q[1];
sx q[1];
rz(1.4742673) q[1];
rz(2.5391606) q[3];
sx q[3];
rz(-0.89655399) q[3];
sx q[3];
rz(-2.7979421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(2.7748761) q[2];
rz(-1.8803053) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(0.096207531) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325608) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(0.60638705) q[0];
rz(2.9442893) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(2.6775449) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7921917) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(2.6130136) q[0];
x q[1];
rz(3.1214141) q[2];
sx q[2];
rz(-1.8115461) q[2];
sx q[2];
rz(-1.1651595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6021767) q[1];
sx q[1];
rz(-2.5578941) q[1];
sx q[1];
rz(1.0023487) q[1];
rz(-pi) q[2];
rz(-0.32902284) q[3];
sx q[3];
rz(-0.24917069) q[3];
sx q[3];
rz(-2.3014625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(-2.8835473) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(-0.095245846) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(1.4415178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1822752) q[0];
sx q[0];
rz(-1.5858985) q[0];
sx q[0];
rz(3.0781151) q[0];
rz(-pi) q[1];
rz(3.0870073) q[2];
sx q[2];
rz(-2.5282801) q[2];
sx q[2];
rz(1.7577946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.085233363) q[1];
sx q[1];
rz(-1.4409522) q[1];
sx q[1];
rz(-1.4443881) q[1];
x q[2];
rz(1.6137245) q[3];
sx q[3];
rz(-1.9044442) q[3];
sx q[3];
rz(1.4294525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(0.22658919) q[2];
rz(-2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-2.3118238) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(-1.7425849) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(-0.98659602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.164455) q[0];
sx q[0];
rz(-1.7286321) q[0];
sx q[0];
rz(-0.49789238) q[0];
rz(-pi) q[1];
rz(-0.95790205) q[2];
sx q[2];
rz(-0.79682486) q[2];
sx q[2];
rz(-2.5755142) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4793195) q[1];
sx q[1];
rz(-0.62367491) q[1];
sx q[1];
rz(-0.24346607) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17296965) q[3];
sx q[3];
rz(-1.9860455) q[3];
sx q[3];
rz(3.1239307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-0.40965664) q[2];
rz(-2.8783197) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.7364527) q[0];
rz(0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.7260889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50354276) q[0];
sx q[0];
rz(-1.8068131) q[0];
sx q[0];
rz(-2.80098) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1591522) q[2];
sx q[2];
rz(-0.42913613) q[2];
sx q[2];
rz(0.8796126) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6092097) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(-1.3551559) q[1];
rz(0.16898445) q[3];
sx q[3];
rz(-2.0770693) q[3];
sx q[3];
rz(2.3224725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(-0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.538095) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(-0.30766906) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(0.60097354) q[2];
sx q[2];
rz(-0.31243159) q[2];
sx q[2];
rz(2.5675788) q[2];
rz(-1.4205167) q[3];
sx q[3];
rz(-1.0192623) q[3];
sx q[3];
rz(-2.7912959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
