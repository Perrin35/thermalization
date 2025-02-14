OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(0.068109186) q[0];
sx q[0];
rz(13.343233) q[0];
rz(-1.3070973) q[1];
sx q[1];
rz(-0.96944648) q[1];
sx q[1];
rz(1.3219272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711219) q[0];
sx q[0];
rz(-1.7092472) q[0];
sx q[0];
rz(-2.2884503) q[0];
rz(-pi) q[1];
rz(0.49643699) q[2];
sx q[2];
rz(-1.8775216) q[2];
sx q[2];
rz(0.53391837) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81404954) q[1];
sx q[1];
rz(-0.40990489) q[1];
sx q[1];
rz(1.4973212) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1573651) q[3];
sx q[3];
rz(-2.7968458) q[3];
sx q[3];
rz(2.1380827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0778568) q[2];
sx q[2];
rz(-1.3602164) q[2];
sx q[2];
rz(2.784101) q[2];
rz(-0.56719559) q[3];
sx q[3];
rz(-1.4128069) q[3];
sx q[3];
rz(-2.7020057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1883063) q[0];
sx q[0];
rz(-0.28306857) q[0];
sx q[0];
rz(-1.8081007) q[0];
rz(-0.4370583) q[1];
sx q[1];
rz(-2.5407365) q[1];
sx q[1];
rz(2.3016047) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8020426) q[0];
sx q[0];
rz(-1.4176482) q[0];
sx q[0];
rz(-2.0846372) q[0];
rz(-pi) q[1];
rz(1.7425547) q[2];
sx q[2];
rz(-2.1585625) q[2];
sx q[2];
rz(-1.5422271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4639791) q[1];
sx q[1];
rz(-2.2871454) q[1];
sx q[1];
rz(0.78732995) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.038957314) q[3];
sx q[3];
rz(-2.1202728) q[3];
sx q[3];
rz(0.4645068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7175032) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(-1.7355512) q[2];
rz(0.16215912) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(-2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9185987) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(-0.42174569) q[0];
rz(-0.84284198) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(0.27580321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3746412) q[0];
sx q[0];
rz(-1.2467442) q[0];
sx q[0];
rz(1.3270017) q[0];
rz(-2.7245443) q[2];
sx q[2];
rz(-2.6337503) q[2];
sx q[2];
rz(-1.5323973) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0581857) q[1];
sx q[1];
rz(-0.32489932) q[1];
sx q[1];
rz(-2.9778381) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86902501) q[3];
sx q[3];
rz(-1.4956253) q[3];
sx q[3];
rz(1.6628671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3383823) q[2];
sx q[2];
rz(-2.3440177) q[2];
sx q[2];
rz(0.69022834) q[2];
rz(2.7817182) q[3];
sx q[3];
rz(-1.3709603) q[3];
sx q[3];
rz(2.7580822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2279219) q[0];
sx q[0];
rz(-1.0424732) q[0];
sx q[0];
rz(0.87164718) q[0];
rz(-0.40920416) q[1];
sx q[1];
rz(-0.49574655) q[1];
sx q[1];
rz(-0.31235487) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.271628) q[0];
sx q[0];
rz(-1.4532491) q[0];
sx q[0];
rz(0.44674504) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4710452) q[2];
sx q[2];
rz(-1.3653697) q[2];
sx q[2];
rz(2.2769986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67955454) q[1];
sx q[1];
rz(-1.7741388) q[1];
sx q[1];
rz(0.36851369) q[1];
x q[2];
rz(1.3757703) q[3];
sx q[3];
rz(-0.93149501) q[3];
sx q[3];
rz(2.2488058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13545869) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(-2.3166166) q[2];
rz(1.7168335) q[3];
sx q[3];
rz(-0.32130876) q[3];
sx q[3];
rz(-2.7062374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.41998) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(-0.87442526) q[0];
rz(-0.32886109) q[1];
sx q[1];
rz(-1.7560274) q[1];
sx q[1];
rz(-1.0985451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.335137) q[0];
sx q[0];
rz(-2.7493434) q[0];
sx q[0];
rz(-0.29229887) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1055036) q[2];
sx q[2];
rz(-0.83613013) q[2];
sx q[2];
rz(-3.0826838) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85531536) q[1];
sx q[1];
rz(-2.554727) q[1];
sx q[1];
rz(2.1333958) q[1];
x q[2];
rz(-2.1608814) q[3];
sx q[3];
rz(-1.3385337) q[3];
sx q[3];
rz(-1.4651608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99981368) q[2];
sx q[2];
rz(-1.1835316) q[2];
sx q[2];
rz(-0.67406526) q[2];
rz(-1.4874124) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(2.8580247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0312626) q[0];
sx q[0];
rz(-2.1188348) q[0];
sx q[0];
rz(1.3856101) q[0];
rz(2.022187) q[1];
sx q[1];
rz(-1.6740572) q[1];
sx q[1];
rz(1.3853692) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35753174) q[0];
sx q[0];
rz(-2.1936532) q[0];
sx q[0];
rz(0.31185598) q[0];
x q[1];
rz(-0.15689705) q[2];
sx q[2];
rz(-0.50373915) q[2];
sx q[2];
rz(2.4258326) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5934903) q[1];
sx q[1];
rz(-1.8917276) q[1];
sx q[1];
rz(1.3671419) q[1];
rz(-pi) q[2];
rz(2.9117111) q[3];
sx q[3];
rz(-2.982989) q[3];
sx q[3];
rz(1.5029217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0783656) q[2];
sx q[2];
rz(-0.93457064) q[2];
sx q[2];
rz(0.37609491) q[2];
rz(-1.1345351) q[3];
sx q[3];
rz(-2.7781656) q[3];
sx q[3];
rz(-2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019808708) q[0];
sx q[0];
rz(-0.60096318) q[0];
sx q[0];
rz(-0.10251775) q[0];
rz(-0.58492297) q[1];
sx q[1];
rz(-2.1939317) q[1];
sx q[1];
rz(1.9580511) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4393504) q[0];
sx q[0];
rz(-0.19225129) q[0];
sx q[0];
rz(-2.4045812) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41924119) q[2];
sx q[2];
rz(-2.4346874) q[2];
sx q[2];
rz(-0.32419328) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.2940604) q[1];
sx q[1];
rz(-2.4449663) q[1];
sx q[1];
rz(-1.9737592) q[1];
x q[2];
rz(-1.3041586) q[3];
sx q[3];
rz(-1.1182866) q[3];
sx q[3];
rz(-3.0686015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3465603) q[2];
sx q[2];
rz(-0.9298032) q[2];
sx q[2];
rz(-2.9417876) q[2];
rz(-2.0810769) q[3];
sx q[3];
rz(-0.54148713) q[3];
sx q[3];
rz(-3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.878433) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(2.8286381) q[0];
rz(-0.9494268) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(1.688028) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87438238) q[0];
sx q[0];
rz(-1.4802762) q[0];
sx q[0];
rz(2.9541914) q[0];
rz(3.0309148) q[2];
sx q[2];
rz(-0.62244906) q[2];
sx q[2];
rz(-1.770293) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18583567) q[1];
sx q[1];
rz(-1.794853) q[1];
sx q[1];
rz(-0.70835374) q[1];
rz(-2.5174934) q[3];
sx q[3];
rz(-0.54275504) q[3];
sx q[3];
rz(-1.8763157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35385418) q[2];
sx q[2];
rz(-0.9757897) q[2];
sx q[2];
rz(-1.806951) q[2];
rz(1.3245964) q[3];
sx q[3];
rz(-2.4748804) q[3];
sx q[3];
rz(-1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1052833) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(-2.4435254) q[0];
rz(2.2659194) q[1];
sx q[1];
rz(-1.061729) q[1];
sx q[1];
rz(0.36044136) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35931906) q[0];
sx q[0];
rz(-2.3397315) q[0];
sx q[0];
rz(2.5129621) q[0];
x q[1];
rz(-2.9827732) q[2];
sx q[2];
rz(-2.7544526) q[2];
sx q[2];
rz(2.5617245) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13294177) q[1];
sx q[1];
rz(-1.1906149) q[1];
sx q[1];
rz(0.93312414) q[1];
rz(-2.8713851) q[3];
sx q[3];
rz(-0.5654618) q[3];
sx q[3];
rz(-2.5266441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91161072) q[2];
sx q[2];
rz(-1.7155827) q[2];
sx q[2];
rz(-2.3243375) q[2];
rz(0.83003712) q[3];
sx q[3];
rz(-2.6006112) q[3];
sx q[3];
rz(0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63477409) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(-3.1095374) q[0];
rz(1.8219148) q[1];
sx q[1];
rz(-1.7664884) q[1];
sx q[1];
rz(2.4874036) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9391032) q[0];
sx q[0];
rz(-1.7669042) q[0];
sx q[0];
rz(1.6926426) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3532969) q[2];
sx q[2];
rz(-1.239778) q[2];
sx q[2];
rz(-2.2904879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75958222) q[1];
sx q[1];
rz(-1.1076265) q[1];
sx q[1];
rz(0.5447556) q[1];
rz(2.8579162) q[3];
sx q[3];
rz(-2.1836238) q[3];
sx q[3];
rz(3.0949288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1303945) q[2];
sx q[2];
rz(-1.0426714) q[2];
sx q[2];
rz(2.1350258) q[2];
rz(-2.448163) q[3];
sx q[3];
rz(-1.3207685) q[3];
sx q[3];
rz(-1.8600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4257767) q[0];
sx q[0];
rz(-2.4632813) q[0];
sx q[0];
rz(3.0671469) q[0];
rz(2.2609932) q[1];
sx q[1];
rz(-1.1136628) q[1];
sx q[1];
rz(-2.640092) q[1];
rz(2.6237189) q[2];
sx q[2];
rz(-2.2020515) q[2];
sx q[2];
rz(-0.9743792) q[2];
rz(2.6941723) q[3];
sx q[3];
rz(-2.3939953) q[3];
sx q[3];
rz(-0.2249997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
