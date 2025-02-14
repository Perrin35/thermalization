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
rz(1.0032049) q[0];
sx q[0];
rz(1.9411074) q[0];
sx q[0];
rz(8.3984126) q[0];
rz(-2.8431471) q[1];
sx q[1];
rz(-1.5084074) q[1];
sx q[1];
rz(2.1389979) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93556773) q[0];
sx q[0];
rz(-1.5479857) q[0];
sx q[0];
rz(-1.3669694) q[0];
rz(-pi) q[1];
rz(0.66291507) q[2];
sx q[2];
rz(-1.6143245) q[2];
sx q[2];
rz(1.7963262) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94384313) q[1];
sx q[1];
rz(-1.437206) q[1];
sx q[1];
rz(0.39033668) q[1];
rz(0.75768023) q[3];
sx q[3];
rz(-2.3509988) q[3];
sx q[3];
rz(2.9942715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58874321) q[2];
sx q[2];
rz(-1.0007977) q[2];
sx q[2];
rz(0.0351077) q[2];
rz(-1.2037753) q[3];
sx q[3];
rz(-1.8193918) q[3];
sx q[3];
rz(-2.8573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0725919) q[0];
sx q[0];
rz(-0.44206107) q[0];
sx q[0];
rz(-2.1373855) q[0];
rz(3.0191811) q[1];
sx q[1];
rz(-1.4675843) q[1];
sx q[1];
rz(0.65707982) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56822424) q[0];
sx q[0];
rz(-1.4889297) q[0];
sx q[0];
rz(0.38827814) q[0];
rz(2.0900656) q[2];
sx q[2];
rz(-2.8247571) q[2];
sx q[2];
rz(0.46665442) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4762558) q[1];
sx q[1];
rz(-1.0968535) q[1];
sx q[1];
rz(-2.302235) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8155926) q[3];
sx q[3];
rz(-2.3395633) q[3];
sx q[3];
rz(2.0662465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.55679) q[2];
sx q[2];
rz(-0.57791296) q[2];
sx q[2];
rz(0.59188262) q[2];
rz(-1.6651734) q[3];
sx q[3];
rz(-1.5341325) q[3];
sx q[3];
rz(0.8736476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0037435) q[0];
sx q[0];
rz(-0.20895222) q[0];
sx q[0];
rz(1.9269706) q[0];
rz(-2.1851723) q[1];
sx q[1];
rz(-1.1914057) q[1];
sx q[1];
rz(-1.1716243) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4223683) q[0];
sx q[0];
rz(-2.5107493) q[0];
sx q[0];
rz(0.8623841) q[0];
rz(-pi) q[1];
rz(-2.9363628) q[2];
sx q[2];
rz(-2.2608888) q[2];
sx q[2];
rz(-0.02366676) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0846935) q[1];
sx q[1];
rz(-2.1325743) q[1];
sx q[1];
rz(1.9719719) q[1];
rz(-pi) q[2];
rz(1.615715) q[3];
sx q[3];
rz(-1.6605596) q[3];
sx q[3];
rz(-2.1819851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9413635) q[2];
sx q[2];
rz(-2.0912632) q[2];
sx q[2];
rz(-1.7639147) q[2];
rz(-2.7348943) q[3];
sx q[3];
rz(-2.4923057) q[3];
sx q[3];
rz(2.6810834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3520626) q[0];
sx q[0];
rz(-1.8822414) q[0];
sx q[0];
rz(-1.9011185) q[0];
rz(0.71818304) q[1];
sx q[1];
rz(-1.3259462) q[1];
sx q[1];
rz(-2.9026418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0485503) q[0];
sx q[0];
rz(-1.6918381) q[0];
sx q[0];
rz(0.28994513) q[0];
x q[1];
rz(-2.4282794) q[2];
sx q[2];
rz(-0.94991604) q[2];
sx q[2];
rz(2.5312405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.093555778) q[1];
sx q[1];
rz(-0.29680064) q[1];
sx q[1];
rz(-1.0637299) q[1];
rz(-pi) q[2];
rz(-0.22438517) q[3];
sx q[3];
rz(-0.54575181) q[3];
sx q[3];
rz(0.70957472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6255528) q[2];
sx q[2];
rz(-1.7480787) q[2];
sx q[2];
rz(-1.3429406) q[2];
rz(-0.88137734) q[3];
sx q[3];
rz(-2.9727029) q[3];
sx q[3];
rz(-2.3605997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4137022) q[0];
sx q[0];
rz(-0.24402937) q[0];
sx q[0];
rz(1.4843041) q[0];
rz(-0.65385747) q[1];
sx q[1];
rz(-1.5212719) q[1];
sx q[1];
rz(-3.003655) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68439517) q[0];
sx q[0];
rz(-0.38551187) q[0];
sx q[0];
rz(1.1657146) q[0];
x q[1];
rz(-2.8213812) q[2];
sx q[2];
rz(-2.4634725) q[2];
sx q[2];
rz(-1.4411639) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0780503) q[1];
sx q[1];
rz(-1.1415328) q[1];
sx q[1];
rz(0.10775609) q[1];
rz(-pi) q[2];
rz(0.055931795) q[3];
sx q[3];
rz(-1.3123296) q[3];
sx q[3];
rz(2.8926769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5518034) q[2];
sx q[2];
rz(-1.9686331) q[2];
sx q[2];
rz(3.0544082) q[2];
rz(0.11463556) q[3];
sx q[3];
rz(-2.3264824) q[3];
sx q[3];
rz(-2.7242928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7755985) q[0];
sx q[0];
rz(-1.5394779) q[0];
sx q[0];
rz(-0.48496801) q[0];
rz(-0.94351774) q[1];
sx q[1];
rz(-2.3285995) q[1];
sx q[1];
rz(-1.2191204) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6935558) q[0];
sx q[0];
rz(-0.72759923) q[0];
sx q[0];
rz(-1.7175402) q[0];
rz(1.2487683) q[2];
sx q[2];
rz(-1.5862984) q[2];
sx q[2];
rz(2.6982234) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.3644219) q[1];
sx q[1];
rz(-0.14816532) q[1];
sx q[1];
rz(1.0309459) q[1];
rz(-pi) q[2];
rz(-0.52690077) q[3];
sx q[3];
rz(-0.8325067) q[3];
sx q[3];
rz(0.47952521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9014827) q[2];
sx q[2];
rz(-2.0811452) q[2];
sx q[2];
rz(0.94129747) q[2];
rz(-2.8955722) q[3];
sx q[3];
rz(-1.4964208) q[3];
sx q[3];
rz(1.7814319) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0435903) q[0];
sx q[0];
rz(-2.0397546) q[0];
sx q[0];
rz(-1.164042) q[0];
rz(-1.3580492) q[1];
sx q[1];
rz(-0.86654228) q[1];
sx q[1];
rz(1.7485626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9214894) q[0];
sx q[0];
rz(-1.1951726) q[0];
sx q[0];
rz(0.079567841) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6955457) q[2];
sx q[2];
rz(-2.5452633) q[2];
sx q[2];
rz(-1.9976384) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23534504) q[1];
sx q[1];
rz(-2.6590912) q[1];
sx q[1];
rz(-0.22884667) q[1];
rz(0.22402899) q[3];
sx q[3];
rz(-2.2371836) q[3];
sx q[3];
rz(-1.9206604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4712269) q[2];
sx q[2];
rz(-1.7134075) q[2];
sx q[2];
rz(-2.8426389) q[2];
rz(0.40545884) q[3];
sx q[3];
rz(-0.36646989) q[3];
sx q[3];
rz(-2.5772742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4295171) q[0];
sx q[0];
rz(-2.8686664) q[0];
sx q[0];
rz(-2.3522229) q[0];
rz(-2.7366267) q[1];
sx q[1];
rz(-1.7908432) q[1];
sx q[1];
rz(0.49496034) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5720737) q[0];
sx q[0];
rz(-1.0099619) q[0];
sx q[0];
rz(0.56604077) q[0];
x q[1];
rz(1.1210006) q[2];
sx q[2];
rz(-0.61857046) q[2];
sx q[2];
rz(0.44396675) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1972754) q[1];
sx q[1];
rz(-2.3305571) q[1];
sx q[1];
rz(1.2326101) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1405419) q[3];
sx q[3];
rz(-2.3214139) q[3];
sx q[3];
rz(-0.4789744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1627545) q[2];
sx q[2];
rz(-1.9163722) q[2];
sx q[2];
rz(-0.2552574) q[2];
rz(-3.1066331) q[3];
sx q[3];
rz(-1.6182263) q[3];
sx q[3];
rz(-2.8953654) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.404945) q[0];
sx q[0];
rz(-1.0628137) q[0];
sx q[0];
rz(2.431562) q[0];
rz(1.0011477) q[1];
sx q[1];
rz(-2.2444057) q[1];
sx q[1];
rz(-0.62090105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9644755) q[0];
sx q[0];
rz(-1.510448) q[0];
sx q[0];
rz(1.2089157) q[0];
rz(-0.59130238) q[2];
sx q[2];
rz(-2.3930051) q[2];
sx q[2];
rz(-0.32574124) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1961851) q[1];
sx q[1];
rz(-1.9176449) q[1];
sx q[1];
rz(1.6849243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.848351) q[3];
sx q[3];
rz(-0.81349868) q[3];
sx q[3];
rz(0.030395776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8093449) q[2];
sx q[2];
rz(-1.5231909) q[2];
sx q[2];
rz(-0.28918949) q[2];
rz(-1.1122164) q[3];
sx q[3];
rz(-2.8465392) q[3];
sx q[3];
rz(2.1212063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4777098) q[0];
sx q[0];
rz(-1.7998671) q[0];
sx q[0];
rz(-2.2148602) q[0];
rz(-1.1070739) q[1];
sx q[1];
rz(-1.5137545) q[1];
sx q[1];
rz(-0.56799299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9290394) q[0];
sx q[0];
rz(-2.4429081) q[0];
sx q[0];
rz(0.50434146) q[0];
rz(-pi) q[1];
x q[1];
rz(1.307042) q[2];
sx q[2];
rz(-2.0677975) q[2];
sx q[2];
rz(-2.618034) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5458249) q[1];
sx q[1];
rz(-2.1819356) q[1];
sx q[1];
rz(0.45054884) q[1];
x q[2];
rz(0.64004489) q[3];
sx q[3];
rz(-1.1403822) q[3];
sx q[3];
rz(1.1120344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0893112) q[2];
sx q[2];
rz(-0.81030446) q[2];
sx q[2];
rz(-1.7029765) q[2];
rz(-0.29664052) q[3];
sx q[3];
rz(-0.34160015) q[3];
sx q[3];
rz(0.44704416) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.457837) q[0];
sx q[0];
rz(-1.0975657) q[0];
sx q[0];
rz(0.56957635) q[0];
rz(2.8029022) q[1];
sx q[1];
rz(-1.6117922) q[1];
sx q[1];
rz(1.502996) q[1];
rz(-1.3647625) q[2];
sx q[2];
rz(-1.3246957) q[2];
sx q[2];
rz(-0.93304721) q[2];
rz(-1.2962616) q[3];
sx q[3];
rz(-1.7332776) q[3];
sx q[3];
rz(0.83040614) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
