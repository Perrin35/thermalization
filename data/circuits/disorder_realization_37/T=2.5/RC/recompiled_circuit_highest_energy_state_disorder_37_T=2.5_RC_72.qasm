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
rz(-1.9395186) q[0];
sx q[0];
rz(2.095686) q[0];
sx q[0];
rz(8.6742824) q[0];
rz(-2.9208288) q[1];
sx q[1];
rz(-1.060744) q[1];
sx q[1];
rz(-1.9579252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25607294) q[0];
sx q[0];
rz(-1.982843) q[0];
sx q[0];
rz(2.4072007) q[0];
rz(-pi) q[1];
rz(-2.2190942) q[2];
sx q[2];
rz(-0.20805173) q[2];
sx q[2];
rz(-1.3626984) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2607076) q[1];
sx q[1];
rz(-1.5591677) q[1];
sx q[1];
rz(2.4616918) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8791844) q[3];
sx q[3];
rz(-1.3779791) q[3];
sx q[3];
rz(-0.28379019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3262647) q[2];
sx q[2];
rz(-1.1417737) q[2];
sx q[2];
rz(-0.092078837) q[2];
rz(-2.0103256) q[3];
sx q[3];
rz(-1.0695894) q[3];
sx q[3];
rz(2.2856975) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6446514) q[0];
sx q[0];
rz(-0.48321378) q[0];
sx q[0];
rz(-0.53571969) q[0];
rz(0.016853111) q[1];
sx q[1];
rz(-2.2622175) q[1];
sx q[1];
rz(0.0171612) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0417174) q[0];
sx q[0];
rz(-1.2493396) q[0];
sx q[0];
rz(-2.0896555) q[0];
x q[1];
rz(-0.13249915) q[2];
sx q[2];
rz(-2.6713058) q[2];
sx q[2];
rz(-2.7072226) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90777724) q[1];
sx q[1];
rz(-1.3388947) q[1];
sx q[1];
rz(1.4091989) q[1];
rz(-3.0073418) q[3];
sx q[3];
rz(-1.6608616) q[3];
sx q[3];
rz(1.7334737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0993593) q[2];
sx q[2];
rz(-1.6681654) q[2];
sx q[2];
rz(0.24822203) q[2];
rz(-0.92784268) q[3];
sx q[3];
rz(-2.35858) q[3];
sx q[3];
rz(2.6718914) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8645653) q[0];
sx q[0];
rz(-1.1772573) q[0];
sx q[0];
rz(-2.9789341) q[0];
rz(-2.3795369) q[1];
sx q[1];
rz(-2.9167852) q[1];
sx q[1];
rz(2.3077097) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57788173) q[0];
sx q[0];
rz(-2.1378008) q[0];
sx q[0];
rz(1.7162697) q[0];
rz(-2.8188848) q[2];
sx q[2];
rz(-1.4526748) q[2];
sx q[2];
rz(-0.28216991) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1385041) q[1];
sx q[1];
rz(-2.2416267) q[1];
sx q[1];
rz(2.6062232) q[1];
x q[2];
rz(-2.7008204) q[3];
sx q[3];
rz(-1.6218054) q[3];
sx q[3];
rz(-0.15404242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0791846) q[2];
sx q[2];
rz(-1.1391613) q[2];
sx q[2];
rz(2.8774234) q[2];
rz(2.0832113) q[3];
sx q[3];
rz(-2.1975785) q[3];
sx q[3];
rz(0.03037608) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77085483) q[0];
sx q[0];
rz(-0.83468947) q[0];
sx q[0];
rz(-1.7133065) q[0];
rz(0.76438534) q[1];
sx q[1];
rz(-1.1420219) q[1];
sx q[1];
rz(-1.3216602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.747731) q[0];
sx q[0];
rz(-1.9698304) q[0];
sx q[0];
rz(1.1927496) q[0];
x q[1];
rz(-1.7197505) q[2];
sx q[2];
rz(-0.33518727) q[2];
sx q[2];
rz(-0.4544979) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.053406) q[1];
sx q[1];
rz(-1.2866396) q[1];
sx q[1];
rz(0.31827815) q[1];
x q[2];
rz(2.0655849) q[3];
sx q[3];
rz(-2.3632999) q[3];
sx q[3];
rz(-0.28631223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4798212) q[2];
sx q[2];
rz(-1.7281374) q[2];
sx q[2];
rz(0.28523764) q[2];
rz(-1.8789004) q[3];
sx q[3];
rz(-1.217239) q[3];
sx q[3];
rz(1.7181905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406463) q[0];
sx q[0];
rz(-0.74746376) q[0];
sx q[0];
rz(-2.2593011) q[0];
rz(-0.67059416) q[1];
sx q[1];
rz(-2.0457485) q[1];
sx q[1];
rz(-0.98463279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4454712) q[0];
sx q[0];
rz(-1.5394256) q[0];
sx q[0];
rz(-0.062792129) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0987344) q[2];
sx q[2];
rz(-1.0295261) q[2];
sx q[2];
rz(1.0554016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.70781103) q[1];
sx q[1];
rz(-1.0278152) q[1];
sx q[1];
rz(2.19853) q[1];
rz(-pi) q[2];
rz(1.0932969) q[3];
sx q[3];
rz(-1.081146) q[3];
sx q[3];
rz(-2.5804105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8661477) q[2];
sx q[2];
rz(-0.53469849) q[2];
sx q[2];
rz(0.61694413) q[2];
rz(-2.3894737) q[3];
sx q[3];
rz(-1.813846) q[3];
sx q[3];
rz(-0.45929685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98944703) q[0];
sx q[0];
rz(-2.4247657) q[0];
sx q[0];
rz(-0.77051198) q[0];
rz(-2.41467) q[1];
sx q[1];
rz(-0.22218552) q[1];
sx q[1];
rz(0.70473421) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4024538) q[0];
sx q[0];
rz(-1.0346892) q[0];
sx q[0];
rz(-2.2741063) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36811916) q[2];
sx q[2];
rz(-2.0621174) q[2];
sx q[2];
rz(-1.4368601) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7721036) q[1];
sx q[1];
rz(-1.9204957) q[1];
sx q[1];
rz(-1.6062801) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9547396) q[3];
sx q[3];
rz(-2.535871) q[3];
sx q[3];
rz(-1.5662409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7695693) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(-2.7083569) q[2];
rz(2.9406934) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(2.6065629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3055426) q[0];
sx q[0];
rz(-1.9925646) q[0];
sx q[0];
rz(-0.38247633) q[0];
rz(-2.1005232) q[1];
sx q[1];
rz(-1.46547) q[1];
sx q[1];
rz(-0.44630757) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6093369) q[0];
sx q[0];
rz(-0.26383886) q[0];
sx q[0];
rz(0.20045073) q[0];
x q[1];
rz(-1.3417753) q[2];
sx q[2];
rz(-1.9441368) q[2];
sx q[2];
rz(2.7315706) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.32532) q[1];
sx q[1];
rz(-2.153646) q[1];
sx q[1];
rz(-2.7144949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71766149) q[3];
sx q[3];
rz(-2.1568686) q[3];
sx q[3];
rz(-1.7505898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3216766) q[2];
sx q[2];
rz(-2.0240462) q[2];
sx q[2];
rz(-0.72810158) q[2];
rz(-0.32650945) q[3];
sx q[3];
rz(-1.6312586) q[3];
sx q[3];
rz(-2.2842893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2364872) q[0];
sx q[0];
rz(-0.40920722) q[0];
sx q[0];
rz(-1.2976728) q[0];
rz(-2.1406651) q[1];
sx q[1];
rz(-0.74556723) q[1];
sx q[1];
rz(2.7330858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7358052) q[0];
sx q[0];
rz(-1.5699495) q[0];
sx q[0];
rz(1.5849831) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3199072) q[2];
sx q[2];
rz(-2.5251881) q[2];
sx q[2];
rz(0.44611713) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7311915) q[1];
sx q[1];
rz(-2.9653314) q[1];
sx q[1];
rz(2.2772863) q[1];
x q[2];
rz(2.767603) q[3];
sx q[3];
rz(-1.9714103) q[3];
sx q[3];
rz(-3.0246322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0268176) q[2];
sx q[2];
rz(-1.6208181) q[2];
sx q[2];
rz(0.82359037) q[2];
rz(-2.1987727) q[3];
sx q[3];
rz(-1.7190245) q[3];
sx q[3];
rz(-1.5028809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5470062) q[0];
sx q[0];
rz(-0.44438812) q[0];
sx q[0];
rz(1.8894926) q[0];
rz(-1.3638672) q[1];
sx q[1];
rz(-1.6412647) q[1];
sx q[1];
rz(1.3346671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9329488) q[0];
sx q[0];
rz(-1.2675537) q[0];
sx q[0];
rz(-1.5327081) q[0];
x q[1];
rz(0.64200114) q[2];
sx q[2];
rz(-1.4307012) q[2];
sx q[2];
rz(-2.7336655) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6393362) q[1];
sx q[1];
rz(-0.4570804) q[1];
sx q[1];
rz(2.9751818) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26504474) q[3];
sx q[3];
rz(-2.0243892) q[3];
sx q[3];
rz(2.3475304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.021726457) q[2];
sx q[2];
rz(-0.29325565) q[2];
sx q[2];
rz(-0.29652706) q[2];
rz(1.4793388) q[3];
sx q[3];
rz(-1.6300423) q[3];
sx q[3];
rz(-0.1571981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5088365) q[0];
sx q[0];
rz(-2.4907676) q[0];
sx q[0];
rz(-0.79488361) q[0];
rz(2.5859313) q[1];
sx q[1];
rz(-1.2763005) q[1];
sx q[1];
rz(1.6937675) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35215615) q[0];
sx q[0];
rz(-1.9481425) q[0];
sx q[0];
rz(0.78363245) q[0];
x q[1];
rz(-3.0270711) q[2];
sx q[2];
rz(-0.50058578) q[2];
sx q[2];
rz(-1.3040744) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8965079) q[1];
sx q[1];
rz(-1.0963529) q[1];
sx q[1];
rz(-2.9175379) q[1];
rz(-pi) q[2];
rz(0.85526222) q[3];
sx q[3];
rz(-2.4088833) q[3];
sx q[3];
rz(-0.39855188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29622233) q[2];
sx q[2];
rz(-1.5105379) q[2];
sx q[2];
rz(1.7108062) q[2];
rz(-0.75302643) q[3];
sx q[3];
rz(-0.81792653) q[3];
sx q[3];
rz(-2.9770765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0059148) q[0];
sx q[0];
rz(-0.68886859) q[0];
sx q[0];
rz(-1.9266358) q[0];
rz(-1.6750492) q[1];
sx q[1];
rz(-0.92862447) q[1];
sx q[1];
rz(-2.5359572) q[1];
rz(0.98943348) q[2];
sx q[2];
rz(-1.0836011) q[2];
sx q[2];
rz(-1.2304162) q[2];
rz(-1.4004272) q[3];
sx q[3];
rz(-1.4650678) q[3];
sx q[3];
rz(2.3847945) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
