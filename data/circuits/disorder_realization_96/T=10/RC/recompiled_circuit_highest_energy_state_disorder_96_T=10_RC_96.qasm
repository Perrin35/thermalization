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
rz(-2.8120256) q[0];
sx q[0];
rz(-0.94533935) q[0];
sx q[0];
rz(3.1402631) q[0];
rz(-2.9381362) q[1];
sx q[1];
rz(-2.9064972) q[1];
sx q[1];
rz(-0.57125616) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093319915) q[0];
sx q[0];
rz(-2.0287345) q[0];
sx q[0];
rz(-2.4620374) q[0];
rz(-2.4680092) q[2];
sx q[2];
rz(-0.24112186) q[2];
sx q[2];
rz(-0.6553638) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23001901) q[1];
sx q[1];
rz(-2.0610448) q[1];
sx q[1];
rz(1.5900299) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9523078) q[3];
sx q[3];
rz(-1.1433868) q[3];
sx q[3];
rz(-1.007146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0433537) q[2];
sx q[2];
rz(-1.2231772) q[2];
sx q[2];
rz(-2.5377048) q[2];
rz(1.0107001) q[3];
sx q[3];
rz(-0.5637919) q[3];
sx q[3];
rz(-0.92196661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0826223) q[0];
sx q[0];
rz(-2.048546) q[0];
sx q[0];
rz(-0.0086722886) q[0];
rz(-1.1588833) q[1];
sx q[1];
rz(-2.4544139) q[1];
sx q[1];
rz(-0.11223665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85189795) q[0];
sx q[0];
rz(-2.0562952) q[0];
sx q[0];
rz(-1.2458891) q[0];
rz(3.0917794) q[2];
sx q[2];
rz(-2.494333) q[2];
sx q[2];
rz(-2.117802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4680926) q[1];
sx q[1];
rz(-1.348494) q[1];
sx q[1];
rz(-2.6373074) q[1];
rz(-pi) q[2];
rz(2.8772164) q[3];
sx q[3];
rz(-1.8473139) q[3];
sx q[3];
rz(-0.087527601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9021641) q[2];
sx q[2];
rz(-1.3210693) q[2];
sx q[2];
rz(-2.2018137) q[2];
rz(2.2043665) q[3];
sx q[3];
rz(-1.7669433) q[3];
sx q[3];
rz(2.7711813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20060191) q[0];
sx q[0];
rz(-0.90018278) q[0];
sx q[0];
rz(-2.5584333) q[0];
rz(1.0519625) q[1];
sx q[1];
rz(-2.5556892) q[1];
sx q[1];
rz(-3.1254056) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28629056) q[0];
sx q[0];
rz(-2.2612766) q[0];
sx q[0];
rz(2.8503522) q[0];
x q[1];
rz(2.3397555) q[2];
sx q[2];
rz(-2.7268598) q[2];
sx q[2];
rz(-1.5842337) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0250108) q[1];
sx q[1];
rz(-1.7336646) q[1];
sx q[1];
rz(-3.0322269) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61516986) q[3];
sx q[3];
rz(-2.2175419) q[3];
sx q[3];
rz(-1.3474479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6215324) q[2];
sx q[2];
rz(-1.3935057) q[2];
sx q[2];
rz(0.072546093) q[2];
rz(2.1324615) q[3];
sx q[3];
rz(-2.3435209) q[3];
sx q[3];
rz(2.9062241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9811454) q[0];
sx q[0];
rz(-0.35780847) q[0];
sx q[0];
rz(-0.88448802) q[0];
rz(1.5011939) q[1];
sx q[1];
rz(-1.6694992) q[1];
sx q[1];
rz(-2.5515058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9035066) q[0];
sx q[0];
rz(-2.0683388) q[0];
sx q[0];
rz(-0.4273129) q[0];
rz(-0.22802148) q[2];
sx q[2];
rz(-1.8393901) q[2];
sx q[2];
rz(-0.10894575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90913016) q[1];
sx q[1];
rz(-1.2386444) q[1];
sx q[1];
rz(-2.0921159) q[1];
rz(1.45118) q[3];
sx q[3];
rz(-1.1362459) q[3];
sx q[3];
rz(-1.0512607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13714743) q[2];
sx q[2];
rz(-1.4819757) q[2];
sx q[2];
rz(1.8972634) q[2];
rz(1.3589877) q[3];
sx q[3];
rz(-2.6748896) q[3];
sx q[3];
rz(0.14314237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2537848) q[0];
sx q[0];
rz(-2.4847327) q[0];
sx q[0];
rz(0.3669056) q[0];
rz(1.3430355) q[1];
sx q[1];
rz(-2.8470706) q[1];
sx q[1];
rz(-1.8273182) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2283233) q[0];
sx q[0];
rz(-2.0073038) q[0];
sx q[0];
rz(0.69605791) q[0];
rz(1.8601244) q[2];
sx q[2];
rz(-2.0380033) q[2];
sx q[2];
rz(-1.6309467) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1340577) q[1];
sx q[1];
rz(-2.4278054) q[1];
sx q[1];
rz(2.0573045) q[1];
x q[2];
rz(2.7845914) q[3];
sx q[3];
rz(-1.9244266) q[3];
sx q[3];
rz(2.2382527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5985976) q[2];
sx q[2];
rz(-2.6948805) q[2];
sx q[2];
rz(-2.561595) q[2];
rz(-0.18481208) q[3];
sx q[3];
rz(-1.5671268) q[3];
sx q[3];
rz(-0.68661657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(3.0762894) q[0];
sx q[0];
rz(-2.6798601) q[0];
sx q[0];
rz(2.8522458) q[0];
rz(0.099960001) q[1];
sx q[1];
rz(-2.3048765) q[1];
sx q[1];
rz(2.3438556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046043175) q[0];
sx q[0];
rz(-2.4975697) q[0];
sx q[0];
rz(-1.796407) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20166534) q[2];
sx q[2];
rz(-1.4603134) q[2];
sx q[2];
rz(0.1486272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5665566) q[1];
sx q[1];
rz(-0.29984353) q[1];
sx q[1];
rz(-2.9672102) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3008505) q[3];
sx q[3];
rz(-1.6463829) q[3];
sx q[3];
rz(-0.10160916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37819698) q[2];
sx q[2];
rz(-1.244647) q[2];
sx q[2];
rz(0.81983105) q[2];
rz(-1.4929006) q[3];
sx q[3];
rz(-0.35455743) q[3];
sx q[3];
rz(-2.7259887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019526871) q[0];
sx q[0];
rz(-2.8391835) q[0];
sx q[0];
rz(2.8269738) q[0];
rz(-2.4954691) q[1];
sx q[1];
rz(-1.4433292) q[1];
sx q[1];
rz(-0.12166469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9537277) q[0];
sx q[0];
rz(-2.6626769) q[0];
sx q[0];
rz(2.2035416) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7971042) q[2];
sx q[2];
rz(-1.0522763) q[2];
sx q[2];
rz(-0.73510494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3683284) q[1];
sx q[1];
rz(-0.46302893) q[1];
sx q[1];
rz(-0.29843389) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8465148) q[3];
sx q[3];
rz(-2.5150617) q[3];
sx q[3];
rz(-2.3008418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6803711) q[2];
sx q[2];
rz(-1.17522) q[2];
sx q[2];
rz(-0.97869527) q[2];
rz(0.77224246) q[3];
sx q[3];
rz(-2.6148836) q[3];
sx q[3];
rz(-2.1286807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1145645) q[0];
sx q[0];
rz(-1.9975198) q[0];
sx q[0];
rz(1.3042599) q[0];
rz(3.0005241) q[1];
sx q[1];
rz(-2.8619659) q[1];
sx q[1];
rz(-1.8581026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1042703) q[0];
sx q[0];
rz(-1.3548343) q[0];
sx q[0];
rz(3.1041935) q[0];
rz(-pi) q[1];
rz(-0.097900585) q[2];
sx q[2];
rz(-1.2719526) q[2];
sx q[2];
rz(3.1063307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5775958) q[1];
sx q[1];
rz(-1.3118) q[1];
sx q[1];
rz(-1.6902655) q[1];
x q[2];
rz(3.0937681) q[3];
sx q[3];
rz(-2.3615814) q[3];
sx q[3];
rz(0.081950233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2804602) q[2];
sx q[2];
rz(-2.2304163) q[2];
sx q[2];
rz(0.99315161) q[2];
rz(-1.6939717) q[3];
sx q[3];
rz(-1.9341035) q[3];
sx q[3];
rz(-1.2710458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3250378) q[0];
sx q[0];
rz(-3.0386381) q[0];
sx q[0];
rz(-0.8763985) q[0];
rz(1.5986298) q[1];
sx q[1];
rz(-1.2878659) q[1];
sx q[1];
rz(2.0230944) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0623613) q[0];
sx q[0];
rz(-0.97506279) q[0];
sx q[0];
rz(-2.0078288) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72965168) q[2];
sx q[2];
rz(-2.0579946) q[2];
sx q[2];
rz(-1.9097569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42414819) q[1];
sx q[1];
rz(-0.70794468) q[1];
sx q[1];
rz(-1.397754) q[1];
rz(-pi) q[2];
rz(0.38409813) q[3];
sx q[3];
rz(-2.4599061) q[3];
sx q[3];
rz(0.1097928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2529651) q[2];
sx q[2];
rz(-0.94433633) q[2];
sx q[2];
rz(-1.2523119) q[2];
rz(-2.0400932) q[3];
sx q[3];
rz(-1.3809729) q[3];
sx q[3];
rz(1.8496877) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0970704) q[0];
sx q[0];
rz(-3.0945393) q[0];
sx q[0];
rz(-0.72730056) q[0];
rz(-2.3749088) q[1];
sx q[1];
rz(-2.2946281) q[1];
sx q[1];
rz(1.1904221) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56404468) q[0];
sx q[0];
rz(-2.2909031) q[0];
sx q[0];
rz(1.4029361) q[0];
x q[1];
rz(1.6929564) q[2];
sx q[2];
rz(-2.4776931) q[2];
sx q[2];
rz(2.4773776) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0990844) q[1];
sx q[1];
rz(-2.7471424) q[1];
sx q[1];
rz(0.69004263) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86807172) q[3];
sx q[3];
rz(-0.70555726) q[3];
sx q[3];
rz(-2.9265273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35655725) q[2];
sx q[2];
rz(-2.0544453) q[2];
sx q[2];
rz(1.1615151) q[2];
rz(-2.0098497) q[3];
sx q[3];
rz(-0.83860207) q[3];
sx q[3];
rz(2.13818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3840735) q[0];
sx q[0];
rz(-2.2390371) q[0];
sx q[0];
rz(-0.83691103) q[0];
rz(1.2548254) q[1];
sx q[1];
rz(-1.246494) q[1];
sx q[1];
rz(1.3424887) q[1];
rz(1.0060476) q[2];
sx q[2];
rz(-1.7972094) q[2];
sx q[2];
rz(-2.1564855) q[2];
rz(0.014160362) q[3];
sx q[3];
rz(-1.7668528) q[3];
sx q[3];
rz(1.7535221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
