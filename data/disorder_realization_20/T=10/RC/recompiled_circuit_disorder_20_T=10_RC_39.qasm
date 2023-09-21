OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.5469172) q[0];
sx q[0];
rz(-2.1407054) q[0];
sx q[0];
rz(0.21240182) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(-1.2815055) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40385383) q[0];
sx q[0];
rz(-1.3699342) q[0];
sx q[0];
rz(1.9181197) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4490657) q[2];
sx q[2];
rz(-1.5480969) q[2];
sx q[2];
rz(-1.9164617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9077397) q[1];
sx q[1];
rz(-2.6820161) q[1];
sx q[1];
rz(1.9563667) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0251861) q[3];
sx q[3];
rz(-2.9985399) q[3];
sx q[3];
rz(1.5798626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4441929) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(-0.15629388) q[0];
rz(-2.5022751) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(1.4888391) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1600906) q[0];
sx q[0];
rz(-2.1182051) q[0];
sx q[0];
rz(-1.7574969) q[0];
rz(-0.82897562) q[2];
sx q[2];
rz(-2.1301221) q[2];
sx q[2];
rz(1.0499357) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7557775) q[1];
sx q[1];
rz(-2.1375244) q[1];
sx q[1];
rz(-2.4789026) q[1];
rz(3.00499) q[3];
sx q[3];
rz(-1.540629) q[3];
sx q[3];
rz(0.68340106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(2.8091649) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333703) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(-2.6112774) q[0];
rz(-2.5976394) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(-1.189032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8408802) q[0];
sx q[0];
rz(-1.8666726) q[0];
sx q[0];
rz(1.4062675) q[0];
x q[1];
rz(-1.5306773) q[2];
sx q[2];
rz(-1.8187858) q[2];
sx q[2];
rz(0.2085533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4424858) q[1];
sx q[1];
rz(-0.92256303) q[1];
sx q[1];
rz(-0.68251619) q[1];
rz(-2.4431908) q[3];
sx q[3];
rz(-1.2050556) q[3];
sx q[3];
rz(1.48667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(0.45423147) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(-2.5355133) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(1.9569424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.143648) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(0.68302897) q[0];
rz(-pi) q[1];
rz(0.84817727) q[2];
sx q[2];
rz(-2.4978673) q[2];
sx q[2];
rz(1.4959469) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6726917) q[1];
sx q[1];
rz(-0.60201445) q[1];
sx q[1];
rz(-1.3637581) q[1];
rz(-2.7167927) q[3];
sx q[3];
rz(-1.9072755) q[3];
sx q[3];
rz(2.6257023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9081395) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(-0.63956368) q[2];
rz(2.284164) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(0.48979428) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63067591) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(-2.6211991) q[0];
rz(3.0386472) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(-2.4127814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6729352) q[0];
sx q[0];
rz(-2.9142083) q[0];
sx q[0];
rz(0.44896455) q[0];
x q[1];
rz(0.16291933) q[2];
sx q[2];
rz(-0.85328057) q[2];
sx q[2];
rz(1.7999072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55620533) q[1];
sx q[1];
rz(-2.7558748) q[1];
sx q[1];
rz(-1.0005887) q[1];
rz(-pi) q[2];
rz(0.95727386) q[3];
sx q[3];
rz(-2.0271218) q[3];
sx q[3];
rz(2.7077655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0985428) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(-0.70971242) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(-3.0098797) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0387886) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(2.2578755) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(-0.19764915) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52299243) q[0];
sx q[0];
rz(-2.1654961) q[0];
sx q[0];
rz(0.80070514) q[0];
rz(-3.0350787) q[2];
sx q[2];
rz(-1.649756) q[2];
sx q[2];
rz(-0.67895141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27053988) q[1];
sx q[1];
rz(-2.6846243) q[1];
sx q[1];
rz(-3.0448826) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88455172) q[3];
sx q[3];
rz(-0.93730799) q[3];
sx q[3];
rz(-1.3214878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(3.0659884) q[2];
rz(1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41912115) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-3.0995195) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(1.1694318) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1202576) q[0];
sx q[0];
rz(-1.406167) q[0];
sx q[0];
rz(3.0740601) q[0];
rz(-1.85497) q[2];
sx q[2];
rz(-1.9467762) q[2];
sx q[2];
rz(0.77501955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3010901) q[1];
sx q[1];
rz(-1.9060578) q[1];
sx q[1];
rz(2.5696978) q[1];
x q[2];
rz(2.7536105) q[3];
sx q[3];
rz(-1.6098607) q[3];
sx q[3];
rz(0.29400533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(0.014483359) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(-0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018205) q[0];
sx q[0];
rz(-0.35035366) q[0];
sx q[0];
rz(-2.5792504) q[0];
rz(1.7100122) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(0.45773488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4938175) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(-1.093822) q[0];
rz(-pi) q[1];
rz(-2.2675603) q[2];
sx q[2];
rz(-2.0119785) q[2];
sx q[2];
rz(0.87768427) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22842912) q[1];
sx q[1];
rz(-2.404681) q[1];
sx q[1];
rz(1.5921028) q[1];
rz(-1.0192972) q[3];
sx q[3];
rz(-2.5428452) q[3];
sx q[3];
rz(2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3056425) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(2.8590554) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(-0.47874513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1829421) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(-1.7392993) q[0];
rz(-2.4422586) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(1.8797849) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14784797) q[0];
sx q[0];
rz(-1.9980668) q[0];
sx q[0];
rz(1.2742313) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4883534) q[2];
sx q[2];
rz(-1.4047033) q[2];
sx q[2];
rz(1.7476029) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3952179) q[1];
sx q[1];
rz(-0.74133855) q[1];
sx q[1];
rz(0.37903255) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7997167) q[3];
sx q[3];
rz(-0.55594) q[3];
sx q[3];
rz(-0.13560175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4385779) q[2];
sx q[2];
rz(-1.6377178) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(2.4750989) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(-0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(0.054140422) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(2.0711526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0265855) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(1.8146145) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1033377) q[2];
sx q[2];
rz(-2.3352211) q[2];
sx q[2];
rz(0.56626608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.479446) q[1];
sx q[1];
rz(-1.3227191) q[1];
sx q[1];
rz(0.76920385) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24786592) q[3];
sx q[3];
rz(-1.7409179) q[3];
sx q[3];
rz(-0.47249139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76876172) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(-1.5268415) q[2];
rz(-0.57957831) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6282745) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(-2.1075481) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(-2.3028159) q[2];
sx q[2];
rz(-1.5427187) q[2];
sx q[2];
rz(-1.6533921) q[2];
rz(-0.53229971) q[3];
sx q[3];
rz(-1.1333864) q[3];
sx q[3];
rz(-0.65896853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];