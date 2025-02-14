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
rz(-2.1276346) q[0];
sx q[0];
rz(-1.4528217) q[0];
sx q[0];
rz(-2.5435574) q[0];
rz(-1.9597837) q[1];
sx q[1];
rz(-2.476517) q[1];
sx q[1];
rz(2.9887078) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476112) q[0];
sx q[0];
rz(-1.6504399) q[0];
sx q[0];
rz(1.682974) q[0];
x q[1];
rz(-3.0316584) q[2];
sx q[2];
rz(-0.49049941) q[2];
sx q[2];
rz(-2.5562499) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8788556) q[1];
sx q[1];
rz(-1.7820184) q[1];
sx q[1];
rz(-3.0131574) q[1];
x q[2];
rz(1.5678207) q[3];
sx q[3];
rz(-1.5248393) q[3];
sx q[3];
rz(-2.9480235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.87024706) q[2];
sx q[2];
rz(-0.18834867) q[2];
sx q[2];
rz(-1.9625473) q[2];
rz(2.7315268) q[3];
sx q[3];
rz(-1.054801) q[3];
sx q[3];
rz(-0.89455354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9533185) q[0];
sx q[0];
rz(-0.66903791) q[0];
sx q[0];
rz(2.8642995) q[0];
rz(-0.20031985) q[1];
sx q[1];
rz(-2.2453313) q[1];
sx q[1];
rz(-3.0576113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6407627) q[0];
sx q[0];
rz(-0.51438722) q[0];
sx q[0];
rz(1.9838833) q[0];
rz(-pi) q[1];
rz(-1.7189156) q[2];
sx q[2];
rz(-0.84858719) q[2];
sx q[2];
rz(2.2715457) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1537495) q[1];
sx q[1];
rz(-1.1293518) q[1];
sx q[1];
rz(2.6325339) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7163058) q[3];
sx q[3];
rz(-0.57568534) q[3];
sx q[3];
rz(3.1120174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1818485) q[2];
sx q[2];
rz(-2.4740348) q[2];
sx q[2];
rz(2.3954771) q[2];
rz(0.6212081) q[3];
sx q[3];
rz(-2.4977903) q[3];
sx q[3];
rz(1.752689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
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
rz(0.32034945) q[0];
sx q[0];
rz(-2.4140883) q[0];
sx q[0];
rz(1.9387091) q[0];
rz(-2.7299643) q[1];
sx q[1];
rz(-1.2806712) q[1];
sx q[1];
rz(-2.0278377) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2611364) q[0];
sx q[0];
rz(-0.87119448) q[0];
sx q[0];
rz(-1.0115795) q[0];
rz(0.97881563) q[2];
sx q[2];
rz(-0.70558682) q[2];
sx q[2];
rz(1.5075114) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74677709) q[1];
sx q[1];
rz(-1.4377497) q[1];
sx q[1];
rz(-0.89292553) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90818543) q[3];
sx q[3];
rz(-2.326283) q[3];
sx q[3];
rz(-0.96409982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0959629) q[2];
sx q[2];
rz(-0.7926422) q[2];
sx q[2];
rz(1.9795798) q[2];
rz(0.80471936) q[3];
sx q[3];
rz(-1.8685124) q[3];
sx q[3];
rz(1.8603676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3941037) q[0];
sx q[0];
rz(-1.2825613) q[0];
sx q[0];
rz(1.5149186) q[0];
rz(2.6331242) q[1];
sx q[1];
rz(-0.6503121) q[1];
sx q[1];
rz(-1.633684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8430458) q[0];
sx q[0];
rz(-0.15104476) q[0];
sx q[0];
rz(-1.5566948) q[0];
rz(-1.3541578) q[2];
sx q[2];
rz(-1.7085953) q[2];
sx q[2];
rz(1.0953643) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.43583664) q[1];
sx q[1];
rz(-1.2661221) q[1];
sx q[1];
rz(3.0235344) q[1];
rz(-pi) q[2];
rz(0.42801492) q[3];
sx q[3];
rz(-2.1819182) q[3];
sx q[3];
rz(1.0304655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7278829) q[2];
sx q[2];
rz(-1.9472313) q[2];
sx q[2];
rz(2.3940562) q[2];
rz(-1.2306635) q[3];
sx q[3];
rz(-0.80941284) q[3];
sx q[3];
rz(0.49720732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.43664765) q[0];
sx q[0];
rz(-1.1282938) q[0];
sx q[0];
rz(-1.4592015) q[0];
rz(-2.7875426) q[1];
sx q[1];
rz(-2.470128) q[1];
sx q[1];
rz(-0.57463247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45978776) q[0];
sx q[0];
rz(-1.5842251) q[0];
sx q[0];
rz(-0.077254967) q[0];
rz(0.76894297) q[2];
sx q[2];
rz(-1.3731602) q[2];
sx q[2];
rz(0.96425024) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3012817) q[1];
sx q[1];
rz(-1.5191829) q[1];
sx q[1];
rz(2.9142996) q[1];
x q[2];
rz(2.0784723) q[3];
sx q[3];
rz(-2.9083544) q[3];
sx q[3];
rz(-0.68768978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6358801) q[2];
sx q[2];
rz(-1.9680223) q[2];
sx q[2];
rz(2.9393348) q[2];
rz(1.0328736) q[3];
sx q[3];
rz(-0.60690108) q[3];
sx q[3];
rz(-0.066298299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.031484) q[0];
sx q[0];
rz(-0.38723543) q[0];
sx q[0];
rz(2.7701344) q[0];
rz(0.29608852) q[1];
sx q[1];
rz(-1.5541872) q[1];
sx q[1];
rz(-0.16337005) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0806769) q[0];
sx q[0];
rz(-1.0759584) q[0];
sx q[0];
rz(-1.7798774) q[0];
rz(-pi) q[1];
rz(-0.54941515) q[2];
sx q[2];
rz(-1.1731469) q[2];
sx q[2];
rz(-2.0009918) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0381546) q[1];
sx q[1];
rz(-1.1327599) q[1];
sx q[1];
rz(-1.7704727) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87518163) q[3];
sx q[3];
rz(-1.4544248) q[3];
sx q[3];
rz(-0.59215183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0801733) q[2];
sx q[2];
rz(-0.27013186) q[2];
sx q[2];
rz(2.488625) q[2];
rz(2.669892) q[3];
sx q[3];
rz(-2.3931849) q[3];
sx q[3];
rz(-2.4145943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93040526) q[0];
sx q[0];
rz(-2.0624332) q[0];
sx q[0];
rz(1.0373254) q[0];
rz(-2.6264181) q[1];
sx q[1];
rz(-1.8637135) q[1];
sx q[1];
rz(-1.3264664) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68149306) q[0];
sx q[0];
rz(-1.1837479) q[0];
sx q[0];
rz(0.99868628) q[0];
rz(-pi) q[1];
rz(0.10978077) q[2];
sx q[2];
rz(-0.35393831) q[2];
sx q[2];
rz(2.372641) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3390914) q[1];
sx q[1];
rz(-0.80858025) q[1];
sx q[1];
rz(1.436961) q[1];
rz(1.8765728) q[3];
sx q[3];
rz(-1.5103087) q[3];
sx q[3];
rz(-0.6636338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6933763) q[2];
sx q[2];
rz(-0.46458149) q[2];
sx q[2];
rz(0.32619897) q[2];
rz(0.97801963) q[3];
sx q[3];
rz(-1.2468485) q[3];
sx q[3];
rz(-3.0514362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4433032) q[0];
sx q[0];
rz(-0.25743085) q[0];
sx q[0];
rz(2.8522016) q[0];
rz(0.012271317) q[1];
sx q[1];
rz(-2.8616276) q[1];
sx q[1];
rz(-1.4853005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0146831) q[0];
sx q[0];
rz(-0.92220798) q[0];
sx q[0];
rz(2.4832583) q[0];
rz(-pi) q[1];
rz(-2.907862) q[2];
sx q[2];
rz(-0.22960358) q[2];
sx q[2];
rz(0.52403852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0694794) q[1];
sx q[1];
rz(-0.24976191) q[1];
sx q[1];
rz(-0.62432557) q[1];
rz(-0.9910219) q[3];
sx q[3];
rz(-2.5977547) q[3];
sx q[3];
rz(2.4286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54899186) q[2];
sx q[2];
rz(-1.5718549) q[2];
sx q[2];
rz(-0.17295095) q[2];
rz(-1.7671827) q[3];
sx q[3];
rz(-1.7632615) q[3];
sx q[3];
rz(1.6636728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7462815) q[0];
sx q[0];
rz(-0.44773856) q[0];
sx q[0];
rz(-2.3315499) q[0];
rz(2.7250302) q[1];
sx q[1];
rz(-2.3655393) q[1];
sx q[1];
rz(-2.7966444) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11703466) q[0];
sx q[0];
rz(-2.2753365) q[0];
sx q[0];
rz(1.7323094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8497963) q[2];
sx q[2];
rz(-1.4115745) q[2];
sx q[2];
rz(2.4375242) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35159416) q[1];
sx q[1];
rz(-1.8934544) q[1];
sx q[1];
rz(2.0568454) q[1];
rz(-pi) q[2];
rz(-0.79699253) q[3];
sx q[3];
rz(-2.9755962) q[3];
sx q[3];
rz(-0.57193631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.714146) q[2];
sx q[2];
rz(-2.2553208) q[2];
sx q[2];
rz(-3.0575338) q[2];
rz(-2.2345624) q[3];
sx q[3];
rz(-1.4182988) q[3];
sx q[3];
rz(0.9575873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8091938) q[0];
sx q[0];
rz(-0.087662307) q[0];
sx q[0];
rz(-2.8293389) q[0];
rz(2.9088083) q[1];
sx q[1];
rz(-2.7504031) q[1];
sx q[1];
rz(-1.2871453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1757619) q[0];
sx q[0];
rz(-1.8928119) q[0];
sx q[0];
rz(-0.67100066) q[0];
rz(-1.5239117) q[2];
sx q[2];
rz(-2.2094036) q[2];
sx q[2];
rz(2.435911) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6617468) q[1];
sx q[1];
rz(-1.3265309) q[1];
sx q[1];
rz(2.0354009) q[1];
x q[2];
rz(2.5165043) q[3];
sx q[3];
rz(-1.2660789) q[3];
sx q[3];
rz(2.8947322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1708258) q[2];
sx q[2];
rz(-1.8584741) q[2];
sx q[2];
rz(-1.4975366) q[2];
rz(2.0434642) q[3];
sx q[3];
rz(-0.69881717) q[3];
sx q[3];
rz(-0.48925492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.004414) q[0];
sx q[0];
rz(-1.351384) q[0];
sx q[0];
rz(2.2030892) q[0];
rz(-1.3068403) q[1];
sx q[1];
rz(-0.88941457) q[1];
sx q[1];
rz(-0.79759146) q[1];
rz(0.14456476) q[2];
sx q[2];
rz(-2.2546386) q[2];
sx q[2];
rz(-2.4972514) q[2];
rz(1.3908006) q[3];
sx q[3];
rz(-0.67435657) q[3];
sx q[3];
rz(-1.6410905) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
