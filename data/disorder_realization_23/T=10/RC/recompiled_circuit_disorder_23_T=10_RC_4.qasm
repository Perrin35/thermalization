OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(-2.5399962) q[1];
sx q[1];
rz(2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9040065) q[0];
sx q[0];
rz(-1.8338086) q[0];
sx q[0];
rz(2.0531274) q[0];
rz(-pi) q[1];
rz(-1.4346052) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(1.9905123) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0093706) q[1];
sx q[1];
rz(-1.8144061) q[1];
sx q[1];
rz(-1.1019215) q[1];
rz(-2.0166287) q[3];
sx q[3];
rz(-0.48005193) q[3];
sx q[3];
rz(-1.1207086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.819954) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(0.83797541) q[2];
rz(-0.49301246) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(-0.078991927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943966) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(-1.863742) q[0];
rz(-0.17678075) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(-0.4321672) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9958187) q[0];
sx q[0];
rz(-0.60129014) q[0];
sx q[0];
rz(-0.50253089) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1439267) q[2];
sx q[2];
rz(-1.9383213) q[2];
sx q[2];
rz(1.2183684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4914815) q[1];
sx q[1];
rz(-1.1307798) q[1];
sx q[1];
rz(-0.12932175) q[1];
rz(-pi) q[2];
rz(3.1167332) q[3];
sx q[3];
rz(-1.0059352) q[3];
sx q[3];
rz(-0.75627518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54923487) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(-0.48669997) q[2];
rz(1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(-2.7242463) q[0];
rz(1.4886645) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(-0.506385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9837512) q[0];
sx q[0];
rz(-2.7451773) q[0];
sx q[0];
rz(1.0979963) q[0];
x q[1];
rz(0.73451368) q[2];
sx q[2];
rz(-1.3909512) q[2];
sx q[2];
rz(-0.70914662) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80458927) q[1];
sx q[1];
rz(-0.5491921) q[1];
sx q[1];
rz(-0.6188436) q[1];
x q[2];
rz(1.9403463) q[3];
sx q[3];
rz(-1.4003716) q[3];
sx q[3];
rz(2.2263262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4425519) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(0.59147269) q[2];
rz(-2.5555723) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(-1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716924) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(-2.3024094) q[0];
rz(0.025578586) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(1.5485839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65028134) q[0];
sx q[0];
rz(-2.6822753) q[0];
sx q[0];
rz(1.7772872) q[0];
rz(-pi) q[1];
rz(2.4315005) q[2];
sx q[2];
rz(-2.2033764) q[2];
sx q[2];
rz(-1.8284423) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7566484) q[1];
sx q[1];
rz(-1.1886016) q[1];
sx q[1];
rz(-0.76809831) q[1];
rz(-pi) q[2];
rz(-1.2876835) q[3];
sx q[3];
rz(-1.1554171) q[3];
sx q[3];
rz(1.1113885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8362391) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(0.099686064) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.8166186) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754159) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(-0.25594041) q[0];
rz(-2.6804965) q[1];
sx q[1];
rz(-2.0979116) q[1];
sx q[1];
rz(-2.3815313) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81239031) q[0];
sx q[0];
rz(-0.15604067) q[0];
sx q[0];
rz(-2.4183256) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7145833) q[2];
sx q[2];
rz(-0.41848768) q[2];
sx q[2];
rz(-0.33868044) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.252994) q[1];
sx q[1];
rz(-2.9610486) q[1];
sx q[1];
rz(0.79043364) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2171668) q[3];
sx q[3];
rz(-1.2741718) q[3];
sx q[3];
rz(0.38277205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(0.6742397) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-3.1242483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6102585) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(-1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(2.0746322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2444006) q[0];
sx q[0];
rz(-0.025248993) q[0];
sx q[0];
rz(-0.61141725) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80337556) q[2];
sx q[2];
rz(-2.5640045) q[2];
sx q[2];
rz(2.9327649) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3031591) q[1];
sx q[1];
rz(-1.618297) q[1];
sx q[1];
rz(-0.82364239) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7339891) q[3];
sx q[3];
rz(-2.3605151) q[3];
sx q[3];
rz(1.2490602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.431488) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(-0.55994326) q[2];
rz(2.4152749) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85754919) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(-1.42111) q[0];
rz(2.9395318) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(0.85817671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7875123) q[0];
sx q[0];
rz(-1.6991827) q[0];
sx q[0];
rz(-1.4228729) q[0];
rz(-pi) q[1];
rz(-3.0965205) q[2];
sx q[2];
rz(-1.9809234) q[2];
sx q[2];
rz(-2.8571667) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8691751) q[1];
sx q[1];
rz(-3.1088447) q[1];
sx q[1];
rz(2.6564024) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5704727) q[3];
sx q[3];
rz(-1.2763599) q[3];
sx q[3];
rz(-1.4049698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28785607) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(1.3254335) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(-0.37471399) q[0];
rz(-0.97887865) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(1.7920866) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5193072) q[0];
sx q[0];
rz(-1.1847727) q[0];
sx q[0];
rz(0.65667721) q[0];
rz(-3.1153203) q[2];
sx q[2];
rz(-0.71825829) q[2];
sx q[2];
rz(-2.9422613) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8229891) q[1];
sx q[1];
rz(-2.9530596) q[1];
sx q[1];
rz(-1.8382501) q[1];
x q[2];
rz(-2.5958063) q[3];
sx q[3];
rz(-1.3457527) q[3];
sx q[3];
rz(3.1080266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65138856) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(-1.7780001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-2.2352495) q[0];
sx q[0];
rz(-1.2619031) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(-1.6040241) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14311929) q[0];
sx q[0];
rz(-2.7612918) q[0];
sx q[0];
rz(-2.2024676) q[0];
rz(0.41861694) q[2];
sx q[2];
rz(-1.4756243) q[2];
sx q[2];
rz(1.0375432) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6919842) q[1];
sx q[1];
rz(-0.87032986) q[1];
sx q[1];
rz(-1.7193754) q[1];
rz(-pi) q[2];
rz(-1.1101301) q[3];
sx q[3];
rz(-1.8970282) q[3];
sx q[3];
rz(-0.5512475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(0.76134479) q[2];
rz(0.90041655) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3392357) q[0];
sx q[0];
rz(-0.46827066) q[0];
sx q[0];
rz(-2.9246869) q[0];
rz(-2.5096109) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(-0.95473081) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.117884) q[0];
sx q[0];
rz(-2.3291991) q[0];
sx q[0];
rz(-0.448416) q[0];
rz(1.9913313) q[2];
sx q[2];
rz(-2.0805801) q[2];
sx q[2];
rz(-1.860294) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4090875) q[1];
sx q[1];
rz(-1.115683) q[1];
sx q[1];
rz(-1.485977) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99777625) q[3];
sx q[3];
rz(-0.7191092) q[3];
sx q[3];
rz(-2.5620808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1251936) q[2];
sx q[2];
rz(-1.23896) q[2];
sx q[2];
rz(-0.94474244) q[2];
rz(-2.7567806) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64086296) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(0.89818556) q[2];
sx q[2];
rz(-2.2834416) q[2];
sx q[2];
rz(-1.7590547) q[2];
rz(2.9363019) q[3];
sx q[3];
rz(-1.1355573) q[3];
sx q[3];
rz(-0.60187403) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
