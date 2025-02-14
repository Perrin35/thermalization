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
rz(0.024700392) q[0];
sx q[0];
rz(-1.2588809) q[0];
sx q[0];
rz(1.2727241) q[0];
rz(1.3920353) q[1];
sx q[1];
rz(-1.4104383) q[1];
sx q[1];
rz(-0.8937723) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0426038) q[0];
sx q[0];
rz(-2.1402142) q[0];
sx q[0];
rz(0.75890394) q[0];
x q[1];
rz(-2.4364901) q[2];
sx q[2];
rz(-0.38422184) q[2];
sx q[2];
rz(-2.9085858) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8302248) q[1];
sx q[1];
rz(-1.1347949) q[1];
sx q[1];
rz(-0.59200759) q[1];
rz(-pi) q[2];
rz(1.5529384) q[3];
sx q[3];
rz(-1.2027338) q[3];
sx q[3];
rz(-2.5758655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36665234) q[2];
sx q[2];
rz(-1.6795936) q[2];
sx q[2];
rz(-2.6095663) q[2];
rz(-0.73412791) q[3];
sx q[3];
rz(-0.2747772) q[3];
sx q[3];
rz(2.0348569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4183913) q[0];
sx q[0];
rz(-2.1558303) q[0];
sx q[0];
rz(-0.080408737) q[0];
rz(-2.6037727) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(0.72371662) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88392576) q[0];
sx q[0];
rz(-2.2307255) q[0];
sx q[0];
rz(1.0038478) q[0];
rz(-2.576648) q[2];
sx q[2];
rz(-2.1568642) q[2];
sx q[2];
rz(2.5069669) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.58077512) q[1];
sx q[1];
rz(-0.58849778) q[1];
sx q[1];
rz(-1.568041) q[1];
rz(-pi) q[2];
rz(-0.81775093) q[3];
sx q[3];
rz(-1.1766542) q[3];
sx q[3];
rz(-0.5241636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85372743) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(2.2354324) q[2];
rz(-0.3624889) q[3];
sx q[3];
rz(-1.5443085) q[3];
sx q[3];
rz(0.046886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.872181) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(2.0500702) q[0];
rz(3.0190234) q[1];
sx q[1];
rz(-1.4361607) q[1];
sx q[1];
rz(1.3465808) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65660921) q[0];
sx q[0];
rz(-0.35728982) q[0];
sx q[0];
rz(1.6553418) q[0];
x q[1];
rz(-2.7287219) q[2];
sx q[2];
rz(-0.70253583) q[2];
sx q[2];
rz(-1.9091878) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31981644) q[1];
sx q[1];
rz(-2.9301661) q[1];
sx q[1];
rz(1.8199304) q[1];
rz(-pi) q[2];
rz(-2.6123206) q[3];
sx q[3];
rz(-1.4733757) q[3];
sx q[3];
rz(-2.4306963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4284105) q[2];
sx q[2];
rz(-1.1298263) q[2];
sx q[2];
rz(-2.909519) q[2];
rz(-1.9706767) q[3];
sx q[3];
rz(-2.5210095) q[3];
sx q[3];
rz(-2.8569729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.93037) q[0];
sx q[0];
rz(-2.3486597) q[0];
sx q[0];
rz(-1.6388182) q[0];
rz(2.5495095) q[1];
sx q[1];
rz(-1.9860257) q[1];
sx q[1];
rz(-0.88964644) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1680206) q[0];
sx q[0];
rz(-2.4102927) q[0];
sx q[0];
rz(-0.84748737) q[0];
rz(-2.0060894) q[2];
sx q[2];
rz(-2.6084628) q[2];
sx q[2];
rz(1.0077623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46737503) q[1];
sx q[1];
rz(-1.2548037) q[1];
sx q[1];
rz(2.9520079) q[1];
x q[2];
rz(2.7569846) q[3];
sx q[3];
rz(-1.9536363) q[3];
sx q[3];
rz(2.2205381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4952116) q[2];
sx q[2];
rz(-1.6782328) q[2];
sx q[2];
rz(-2.4656673) q[2];
rz(1.4365139) q[3];
sx q[3];
rz(-0.33994514) q[3];
sx q[3];
rz(-2.9743312) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2109569) q[0];
sx q[0];
rz(-2.8093331) q[0];
sx q[0];
rz(-0.19530547) q[0];
rz(-3.0209814) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(-2.9511071) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2094177) q[0];
sx q[0];
rz(-0.10917347) q[0];
sx q[0];
rz(2.5072844) q[0];
rz(-pi) q[1];
rz(1.6570857) q[2];
sx q[2];
rz(-0.89601529) q[2];
sx q[2];
rz(2.9099885) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45791679) q[1];
sx q[1];
rz(-1.6714339) q[1];
sx q[1];
rz(1.5186054) q[1];
rz(-pi) q[2];
rz(0.61483947) q[3];
sx q[3];
rz(-1.9214464) q[3];
sx q[3];
rz(-0.21803724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5087937) q[2];
sx q[2];
rz(-1.4845279) q[2];
sx q[2];
rz(-1.9602027) q[2];
rz(1.7975636) q[3];
sx q[3];
rz(-2.400178) q[3];
sx q[3];
rz(-0.77176362) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33264273) q[0];
sx q[0];
rz(-1.7985666) q[0];
sx q[0];
rz(-1.1395662) q[0];
rz(-0.66967669) q[1];
sx q[1];
rz(-2.2256336) q[1];
sx q[1];
rz(-2.0119827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87605873) q[0];
sx q[0];
rz(-1.8647412) q[0];
sx q[0];
rz(-2.8493899) q[0];
rz(-pi) q[1];
rz(-1.3048346) q[2];
sx q[2];
rz(-2.2707191) q[2];
sx q[2];
rz(-1.0476607) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5796389) q[1];
sx q[1];
rz(-0.41746751) q[1];
sx q[1];
rz(-2.6446656) q[1];
x q[2];
rz(2.4242006) q[3];
sx q[3];
rz(-2.903707) q[3];
sx q[3];
rz(-0.57008509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.42664042) q[2];
sx q[2];
rz(-1.7317829) q[2];
sx q[2];
rz(-0.40531522) q[2];
rz(-2.3146368) q[3];
sx q[3];
rz(-0.6260286) q[3];
sx q[3];
rz(2.2820182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7760794) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(-2.205701) q[0];
rz(0.86839688) q[1];
sx q[1];
rz(-2.6135542) q[1];
sx q[1];
rz(-1.7220928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065212981) q[0];
sx q[0];
rz(-1.924949) q[0];
sx q[0];
rz(0.54535086) q[0];
x q[1];
rz(-2.9090857) q[2];
sx q[2];
rz(-1.075313) q[2];
sx q[2];
rz(-1.8956309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81379997) q[1];
sx q[1];
rz(-2.2918502) q[1];
sx q[1];
rz(2.2071597) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0017599864) q[3];
sx q[3];
rz(-1.5704182) q[3];
sx q[3];
rz(0.019824713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.821227) q[2];
sx q[2];
rz(-1.3996404) q[2];
sx q[2];
rz(2.121675) q[2];
rz(2.6323281) q[3];
sx q[3];
rz(-1.7002707) q[3];
sx q[3];
rz(0.078484623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34198636) q[0];
sx q[0];
rz(-1.5155563) q[0];
sx q[0];
rz(-1.9827783) q[0];
rz(-0.47053567) q[1];
sx q[1];
rz(-1.7262986) q[1];
sx q[1];
rz(1.4656969) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19150951) q[0];
sx q[0];
rz(-1.5198764) q[0];
sx q[0];
rz(-1.5470734) q[0];
rz(2.450804) q[2];
sx q[2];
rz(-1.2993408) q[2];
sx q[2];
rz(-1.1784306) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.012124) q[1];
sx q[1];
rz(-0.91255847) q[1];
sx q[1];
rz(1.6679428) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37181446) q[3];
sx q[3];
rz(-0.68228693) q[3];
sx q[3];
rz(2.2390847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1450242) q[2];
sx q[2];
rz(-1.8137167) q[2];
sx q[2];
rz(-2.3222951) q[2];
rz(0.79646349) q[3];
sx q[3];
rz(-2.7345246) q[3];
sx q[3];
rz(1.4888633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2040937) q[0];
sx q[0];
rz(-0.096385328) q[0];
sx q[0];
rz(2.3199484) q[0];
rz(2.4687528) q[1];
sx q[1];
rz(-0.57603374) q[1];
sx q[1];
rz(-1.0427262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0053528) q[0];
sx q[0];
rz(-1.3658821) q[0];
sx q[0];
rz(-1.4968275) q[0];
x q[1];
rz(-0.64260245) q[2];
sx q[2];
rz(-2.2589141) q[2];
sx q[2];
rz(-2.7049899) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71127001) q[1];
sx q[1];
rz(-1.9414328) q[1];
sx q[1];
rz(2.1881585) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.72307) q[3];
sx q[3];
rz(-1.8209753) q[3];
sx q[3];
rz(-2.6891249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3599856) q[2];
sx q[2];
rz(-1.8693417) q[2];
sx q[2];
rz(0.84235111) q[2];
rz(-0.53884566) q[3];
sx q[3];
rz(-1.6539961) q[3];
sx q[3];
rz(-2.6174788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46776029) q[0];
sx q[0];
rz(-0.14586511) q[0];
sx q[0];
rz(1.9061331) q[0];
rz(2.1927059) q[1];
sx q[1];
rz(-2.3088539) q[1];
sx q[1];
rz(-1.4804776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6687889) q[0];
sx q[0];
rz(-1.2929217) q[0];
sx q[0];
rz(-0.86937277) q[0];
rz(-pi) q[1];
rz(-0.95959227) q[2];
sx q[2];
rz(-1.094154) q[2];
sx q[2];
rz(0.075919064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9292641) q[1];
sx q[1];
rz(-2.9605299) q[1];
sx q[1];
rz(-2.7140815) q[1];
x q[2];
rz(0.62678316) q[3];
sx q[3];
rz(-2.7434182) q[3];
sx q[3];
rz(-2.6836723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4296055) q[2];
sx q[2];
rz(-1.4236071) q[2];
sx q[2];
rz(-2.831366) q[2];
rz(2.6817536) q[3];
sx q[3];
rz(-0.55791563) q[3];
sx q[3];
rz(1.7673813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9643758) q[0];
sx q[0];
rz(-2.0073267) q[0];
sx q[0];
rz(1.0649756) q[0];
rz(-2.2433157) q[1];
sx q[1];
rz(-2.5986462) q[1];
sx q[1];
rz(2.6959261) q[1];
rz(-0.40681268) q[2];
sx q[2];
rz(-2.1944703) q[2];
sx q[2];
rz(-2.2866972) q[2];
rz(2.0528005) q[3];
sx q[3];
rz(-0.72127753) q[3];
sx q[3];
rz(0.75281561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
