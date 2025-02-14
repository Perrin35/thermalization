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
rz(0.32956707) q[0];
sx q[0];
rz(0.94533935) q[0];
sx q[0];
rz(9.4234484) q[0];
rz(-2.9381362) q[1];
sx q[1];
rz(-2.9064972) q[1];
sx q[1];
rz(-0.57125616) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0482727) q[0];
sx q[0];
rz(-2.0287345) q[0];
sx q[0];
rz(0.67955525) q[0];
rz(2.9517101) q[2];
sx q[2];
rz(-1.7203091) q[2];
sx q[2];
rz(0.25611311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9115736) q[1];
sx q[1];
rz(-2.0610448) q[1];
sx q[1];
rz(1.5900299) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6853525) q[3];
sx q[3];
rz(-1.9164698) q[3];
sx q[3];
rz(0.72843188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0433537) q[2];
sx q[2];
rz(-1.9184155) q[2];
sx q[2];
rz(2.5377048) q[2];
rz(-2.1308925) q[3];
sx q[3];
rz(-2.5778008) q[3];
sx q[3];
rz(-2.219626) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0589703) q[0];
sx q[0];
rz(-2.048546) q[0];
sx q[0];
rz(-0.0086722886) q[0];
rz(-1.9827093) q[1];
sx q[1];
rz(-2.4544139) q[1];
sx q[1];
rz(0.11223665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4771292) q[0];
sx q[0];
rz(-2.5647129) q[0];
sx q[0];
rz(0.5441027) q[0];
x q[1];
rz(2.49493) q[2];
sx q[2];
rz(-1.5407667) q[2];
sx q[2];
rz(0.5867556) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.018377233) q[1];
sx q[1];
rz(-2.0615512) q[1];
sx q[1];
rz(-1.823455) q[1];
rz(-1.2848507) q[3];
sx q[3];
rz(-1.8249028) q[3];
sx q[3];
rz(1.4094963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.23942854) q[2];
sx q[2];
rz(-1.8205234) q[2];
sx q[2];
rz(0.93977896) q[2];
rz(-0.93722614) q[3];
sx q[3];
rz(-1.3746494) q[3];
sx q[3];
rz(-2.7711813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9409907) q[0];
sx q[0];
rz(-2.2414099) q[0];
sx q[0];
rz(-2.5584333) q[0];
rz(-2.0896301) q[1];
sx q[1];
rz(-2.5556892) q[1];
sx q[1];
rz(0.016187035) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72618851) q[0];
sx q[0];
rz(-2.4016018) q[0];
sx q[0];
rz(1.9053024) q[0];
rz(-0.80183713) q[2];
sx q[2];
rz(-0.41473284) q[2];
sx q[2];
rz(1.5842337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5213373) q[1];
sx q[1];
rz(-2.9456821) q[1];
sx q[1];
rz(-2.1571062) q[1];
x q[2];
rz(-0.61516986) q[3];
sx q[3];
rz(-0.92405073) q[3];
sx q[3];
rz(-1.3474479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6215324) q[2];
sx q[2];
rz(-1.3935057) q[2];
sx q[2];
rz(-0.072546093) q[2];
rz(-1.0091311) q[3];
sx q[3];
rz(-0.79807177) q[3];
sx q[3];
rz(0.23536853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16044727) q[0];
sx q[0];
rz(-2.7837842) q[0];
sx q[0];
rz(2.2571046) q[0];
rz(-1.6403987) q[1];
sx q[1];
rz(-1.6694992) q[1];
sx q[1];
rz(-2.5515058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4761596) q[0];
sx q[0];
rz(-2.4976625) q[0];
sx q[0];
rz(2.2225999) q[0];
rz(-pi) q[1];
rz(2.2583738) q[2];
sx q[2];
rz(-0.35053634) q[2];
sx q[2];
rz(0.8274629) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47653884) q[1];
sx q[1];
rz(-2.0609984) q[1];
sx q[1];
rz(-2.7630105) q[1];
x q[2];
rz(-2.7042974) q[3];
sx q[3];
rz(-1.4623433) q[3];
sx q[3];
rz(-2.5714998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13714743) q[2];
sx q[2];
rz(-1.4819757) q[2];
sx q[2];
rz(1.8972634) q[2];
rz(-1.7826049) q[3];
sx q[3];
rz(-2.6748896) q[3];
sx q[3];
rz(0.14314237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8878079) q[0];
sx q[0];
rz(-2.4847327) q[0];
sx q[0];
rz(0.3669056) q[0];
rz(1.7985571) q[1];
sx q[1];
rz(-2.8470706) q[1];
sx q[1];
rz(-1.3142745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99710354) q[0];
sx q[0];
rz(-2.1908452) q[0];
sx q[0];
rz(-1.0245566) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6268705) q[2];
sx q[2];
rz(-2.597737) q[2];
sx q[2];
rz(-1.0469231) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.007535) q[1];
sx q[1];
rz(-2.4278054) q[1];
sx q[1];
rz(-1.0842881) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7845914) q[3];
sx q[3];
rz(-1.9244266) q[3];
sx q[3];
rz(2.2382527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.542995) q[2];
sx q[2];
rz(-2.6948805) q[2];
sx q[2];
rz(0.57999769) q[2];
rz(-0.18481208) q[3];
sx q[3];
rz(-1.5744659) q[3];
sx q[3];
rz(0.68661657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065303236) q[0];
sx q[0];
rz(-0.46173254) q[0];
sx q[0];
rz(2.8522458) q[0];
rz(-0.099960001) q[1];
sx q[1];
rz(-0.83671612) q[1];
sx q[1];
rz(-0.79773703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3432309) q[0];
sx q[0];
rz(-1.7055178) q[0];
sx q[0];
rz(0.93905296) q[0];
rz(-pi) q[1];
rz(1.6835454) q[2];
sx q[2];
rz(-1.7712153) q[2];
sx q[2];
rz(-1.7419614) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7489075) q[1];
sx q[1];
rz(-1.2756383) q[1];
sx q[1];
rz(-1.6243851) q[1];
x q[2];
rz(-2.3008505) q[3];
sx q[3];
rz(-1.6463829) q[3];
sx q[3];
rz(0.10160916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7633957) q[2];
sx q[2];
rz(-1.8969456) q[2];
sx q[2];
rz(2.3217616) q[2];
rz(-1.648692) q[3];
sx q[3];
rz(-2.7870352) q[3];
sx q[3];
rz(-2.7259887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019526871) q[0];
sx q[0];
rz(-2.8391835) q[0];
sx q[0];
rz(-0.31461883) q[0];
rz(-0.64612359) q[1];
sx q[1];
rz(-1.4433292) q[1];
sx q[1];
rz(0.12166469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3356161) q[0];
sx q[0];
rz(-1.2947963) q[0];
sx q[0];
rz(1.9673303) q[0];
rz(-pi) q[1];
rz(-0.7971042) q[2];
sx q[2];
rz(-1.0522763) q[2];
sx q[2];
rz(2.4064877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52895228) q[1];
sx q[1];
rz(-1.4390872) q[1];
sx q[1];
rz(-0.44520933) q[1];
rz(-pi) q[2];
rz(-2.9470351) q[3];
sx q[3];
rz(-0.97134198) q[3];
sx q[3];
rz(-0.50473467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6803711) q[2];
sx q[2];
rz(-1.9663726) q[2];
sx q[2];
rz(0.97869527) q[2];
rz(-0.77224246) q[3];
sx q[3];
rz(-2.6148836) q[3];
sx q[3];
rz(-1.0129119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-1.1145645) q[0];
sx q[0];
rz(-1.1440729) q[0];
sx q[0];
rz(1.3042599) q[0];
rz(-0.14106855) q[1];
sx q[1];
rz(-2.8619659) q[1];
sx q[1];
rz(1.28349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9314037) q[0];
sx q[0];
rz(-0.21912665) q[0];
sx q[0];
rz(1.7396084) q[0];
x q[1];
rz(-1.2705994) q[2];
sx q[2];
rz(-1.4772479) q[2];
sx q[2];
rz(1.5644422) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.12569553) q[1];
sx q[1];
rz(-0.28465965) q[1];
sx q[1];
rz(2.7188676) q[1];
rz(-2.3621534) q[3];
sx q[3];
rz(-1.6044242) q[3];
sx q[3];
rz(1.5228576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2804602) q[2];
sx q[2];
rz(-2.2304163) q[2];
sx q[2];
rz(0.99315161) q[2];
rz(-1.447621) q[3];
sx q[3];
rz(-1.2074892) q[3];
sx q[3];
rz(1.8705468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3250378) q[0];
sx q[0];
rz(-0.10295454) q[0];
sx q[0];
rz(2.2651941) q[0];
rz(1.5986298) q[1];
sx q[1];
rz(-1.8537268) q[1];
sx q[1];
rz(1.1184982) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23520175) q[0];
sx q[0];
rz(-1.9287325) q[0];
sx q[0];
rz(-2.4992348) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72965168) q[2];
sx q[2];
rz(-1.083598) q[2];
sx q[2];
rz(-1.9097569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4913018) q[1];
sx q[1];
rz(-0.87556616) q[1];
sx q[1];
rz(0.14632605) q[1];
rz(-pi) q[2];
rz(-0.64499027) q[3];
sx q[3];
rz(-1.8091615) q[3];
sx q[3];
rz(1.1568943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8886275) q[2];
sx q[2];
rz(-0.94433633) q[2];
sx q[2];
rz(-1.8892807) q[2];
rz(2.0400932) q[3];
sx q[3];
rz(-1.7606198) q[3];
sx q[3];
rz(-1.291905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0445223) q[0];
sx q[0];
rz(-0.047053311) q[0];
sx q[0];
rz(0.72730056) q[0];
rz(-2.3749088) q[1];
sx q[1];
rz(-2.2946281) q[1];
sx q[1];
rz(-1.9511706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81556126) q[0];
sx q[0];
rz(-2.4056068) q[0];
sx q[0];
rz(2.953397) q[0];
x q[1];
rz(-1.6929564) q[2];
sx q[2];
rz(-2.4776931) q[2];
sx q[2];
rz(-2.4773776) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3695581) q[1];
sx q[1];
rz(-1.2698962) q[1];
sx q[1];
rz(1.3117666) q[1];
rz(-pi) q[2];
rz(-2.2735209) q[3];
sx q[3];
rz(-0.70555726) q[3];
sx q[3];
rz(-0.21506532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35655725) q[2];
sx q[2];
rz(-1.0871474) q[2];
sx q[2];
rz(-1.9800775) q[2];
rz(2.0098497) q[3];
sx q[3];
rz(-2.3029906) q[3];
sx q[3];
rz(-1.0034126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.9772498) q[2];
sx q[2];
rz(-0.60383527) q[2];
sx q[2];
rz(2.8964103) q[2];
rz(-3.1274323) q[3];
sx q[3];
rz(-1.7668528) q[3];
sx q[3];
rz(1.7535221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
