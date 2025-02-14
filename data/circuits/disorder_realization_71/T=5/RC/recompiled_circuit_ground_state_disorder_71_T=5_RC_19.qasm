OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.355964) q[0];
sx q[0];
rz(-2.845919) q[0];
sx q[0];
rz(-0.41391882) q[0];
rz(3.4105372) q[1];
sx q[1];
rz(6.6528448) q[1];
sx q[1];
rz(11.289968) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9577173) q[0];
sx q[0];
rz(-1.661429) q[0];
sx q[0];
rz(3.0498226) q[0];
rz(1.4914614) q[2];
sx q[2];
rz(-0.66986194) q[2];
sx q[2];
rz(1.8928493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1943097) q[1];
sx q[1];
rz(-2.3230246) q[1];
sx q[1];
rz(1.964393) q[1];
x q[2];
rz(1.3876565) q[3];
sx q[3];
rz(-2.1570342) q[3];
sx q[3];
rz(-0.30979118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3966177) q[2];
sx q[2];
rz(-1.2367542) q[2];
sx q[2];
rz(-1.940894) q[2];
rz(-2.4602304) q[3];
sx q[3];
rz(-1.6187637) q[3];
sx q[3];
rz(2.7549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6354527) q[0];
sx q[0];
rz(-2.8400087) q[0];
sx q[0];
rz(-0.26148456) q[0];
rz(1.5226978) q[1];
sx q[1];
rz(-0.50917429) q[1];
sx q[1];
rz(1.8656628) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44156333) q[0];
sx q[0];
rz(-1.4744548) q[0];
sx q[0];
rz(-0.55192134) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4894756) q[2];
sx q[2];
rz(-2.0154833) q[2];
sx q[2];
rz(-2.7114078) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51569444) q[1];
sx q[1];
rz(-1.0401023) q[1];
sx q[1];
rz(-1.9489392) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1072282) q[3];
sx q[3];
rz(-1.4219163) q[3];
sx q[3];
rz(1.9158165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2168938) q[2];
sx q[2];
rz(-1.0328707) q[2];
sx q[2];
rz(0.72803289) q[2];
rz(-1.3701471) q[3];
sx q[3];
rz(-0.61714211) q[3];
sx q[3];
rz(-0.035577687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034123357) q[0];
sx q[0];
rz(-1.1119482) q[0];
sx q[0];
rz(0.6849826) q[0];
rz(-2.0383539) q[1];
sx q[1];
rz(-0.61956844) q[1];
sx q[1];
rz(-1.0122976) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2489126) q[0];
sx q[0];
rz(-1.9143701) q[0];
sx q[0];
rz(-1.1569886) q[0];
rz(-pi) q[1];
rz(2.9673789) q[2];
sx q[2];
rz(-1.763134) q[2];
sx q[2];
rz(-1.9680374) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5372181) q[1];
sx q[1];
rz(-2.8979725) q[1];
sx q[1];
rz(2.591406) q[1];
x q[2];
rz(2.3776953) q[3];
sx q[3];
rz(-1.9490846) q[3];
sx q[3];
rz(2.0495142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.181902) q[2];
sx q[2];
rz(-1.5519374) q[2];
sx q[2];
rz(2.336179) q[2];
rz(2.3253333) q[3];
sx q[3];
rz(-0.62556848) q[3];
sx q[3];
rz(0.094154112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020298088) q[0];
sx q[0];
rz(-0.90885201) q[0];
sx q[0];
rz(0.62776172) q[0];
rz(-3.032387) q[1];
sx q[1];
rz(-2.2716378) q[1];
sx q[1];
rz(0.83762082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254346) q[0];
sx q[0];
rz(-1.6372879) q[0];
sx q[0];
rz(2.6496135) q[0];
x q[1];
rz(-1.8960647) q[2];
sx q[2];
rz(-0.62042371) q[2];
sx q[2];
rz(1.9577946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.38214907) q[1];
sx q[1];
rz(-1.4370222) q[1];
sx q[1];
rz(0.70258883) q[1];
rz(2.4627962) q[3];
sx q[3];
rz(-0.46884109) q[3];
sx q[3];
rz(-1.9285338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8236905) q[2];
sx q[2];
rz(-2.6573942) q[2];
sx q[2];
rz(-2.79706) q[2];
rz(-1.4065929) q[3];
sx q[3];
rz(-0.35511261) q[3];
sx q[3];
rz(0.047903456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.542881) q[0];
sx q[0];
rz(-0.6364091) q[0];
sx q[0];
rz(-2.6357546) q[0];
rz(-2.089031) q[1];
sx q[1];
rz(-2.274175) q[1];
sx q[1];
rz(2.1956086) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3179683) q[0];
sx q[0];
rz(-1.8232947) q[0];
sx q[0];
rz(0.28244762) q[0];
rz(-1.4267816) q[2];
sx q[2];
rz(-1.9249232) q[2];
sx q[2];
rz(1.9812552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8648423) q[1];
sx q[1];
rz(-0.71886834) q[1];
sx q[1];
rz(2.5548773) q[1];
rz(-0.068377062) q[3];
sx q[3];
rz(-1.42236) q[3];
sx q[3];
rz(-2.1297586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0537009) q[2];
sx q[2];
rz(-1.0975857) q[2];
sx q[2];
rz(-1.5778479) q[2];
rz(2.0569233) q[3];
sx q[3];
rz(-0.8756777) q[3];
sx q[3];
rz(-2.5904371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56187335) q[0];
sx q[0];
rz(-2.5358574) q[0];
sx q[0];
rz(-2.9015923) q[0];
rz(-2.2198246) q[1];
sx q[1];
rz(-2.0926937) q[1];
sx q[1];
rz(1.1704495) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76239785) q[0];
sx q[0];
rz(-1.0012915) q[0];
sx q[0];
rz(-1.5194917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7851849) q[2];
sx q[2];
rz(-1.6898167) q[2];
sx q[2];
rz(-2.1343892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6854418) q[1];
sx q[1];
rz(-2.4439881) q[1];
sx q[1];
rz(-0.82932034) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0441566) q[3];
sx q[3];
rz(-1.9698922) q[3];
sx q[3];
rz(-2.4533085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74448284) q[2];
sx q[2];
rz(-0.8195256) q[2];
sx q[2];
rz(-2.5605555) q[2];
rz(0.93823141) q[3];
sx q[3];
rz(-2.5751028) q[3];
sx q[3];
rz(-0.70403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0051603) q[0];
sx q[0];
rz(-0.67661023) q[0];
sx q[0];
rz(-2.6410979) q[0];
rz(-0.23807921) q[1];
sx q[1];
rz(-1.6903189) q[1];
sx q[1];
rz(-0.64661017) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4256) q[0];
sx q[0];
rz(-1.5725721) q[0];
sx q[0];
rz(3.1067757) q[0];
rz(3.1076381) q[2];
sx q[2];
rz(-1.6852549) q[2];
sx q[2];
rz(-1.7781228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2518066) q[1];
sx q[1];
rz(-1.6908892) q[1];
sx q[1];
rz(-1.0838064) q[1];
rz(-1.2122596) q[3];
sx q[3];
rz(-1.6286116) q[3];
sx q[3];
rz(2.0023579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2639192) q[2];
sx q[2];
rz(-0.71920243) q[2];
sx q[2];
rz(-0.70362299) q[2];
rz(-1.8536812) q[3];
sx q[3];
rz(-1.0602919) q[3];
sx q[3];
rz(-2.6357486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8877761) q[0];
sx q[0];
rz(-1.1340589) q[0];
sx q[0];
rz(-1.3251086) q[0];
rz(-2.3690986) q[1];
sx q[1];
rz(-1.1224116) q[1];
sx q[1];
rz(-2.966029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2529375) q[0];
sx q[0];
rz(-2.3217259) q[0];
sx q[0];
rz(1.8273414) q[0];
rz(-2.3374723) q[2];
sx q[2];
rz(-2.0607227) q[2];
sx q[2];
rz(2.8049289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0766636) q[1];
sx q[1];
rz(-2.6573219) q[1];
sx q[1];
rz(-1.2679497) q[1];
rz(2.9977071) q[3];
sx q[3];
rz(-2.4200614) q[3];
sx q[3];
rz(0.084621457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10302155) q[2];
sx q[2];
rz(-2.0671637) q[2];
sx q[2];
rz(1.0938905) q[2];
rz(-1.966018) q[3];
sx q[3];
rz(-0.75785494) q[3];
sx q[3];
rz(-0.28483835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71838251) q[0];
sx q[0];
rz(-0.045561401) q[0];
sx q[0];
rz(1.3257931) q[0];
rz(-2.6153053) q[1];
sx q[1];
rz(-0.86295366) q[1];
sx q[1];
rz(-2.3525499) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51368749) q[0];
sx q[0];
rz(-2.9531678) q[0];
sx q[0];
rz(1.0303251) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7858235) q[2];
sx q[2];
rz(-1.8104743) q[2];
sx q[2];
rz(0.046047839) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0328417) q[1];
sx q[1];
rz(-0.95966731) q[1];
sx q[1];
rz(-2.717589) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2317722) q[3];
sx q[3];
rz(-2.2299521) q[3];
sx q[3];
rz(0.59162635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9290756) q[2];
sx q[2];
rz(-0.93364659) q[2];
sx q[2];
rz(2.2361501) q[2];
rz(2.2255911) q[3];
sx q[3];
rz(-1.4429049) q[3];
sx q[3];
rz(2.0974832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4395897) q[0];
sx q[0];
rz(-2.3083394) q[0];
sx q[0];
rz(-0.92779094) q[0];
rz(1.0935498) q[1];
sx q[1];
rz(-0.55892006) q[1];
sx q[1];
rz(0.13339001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0463842) q[0];
sx q[0];
rz(-2.781087) q[0];
sx q[0];
rz(2.5166839) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8850408) q[2];
sx q[2];
rz(-2.5470599) q[2];
sx q[2];
rz(2.9192215) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3747537) q[1];
sx q[1];
rz(-2.1275318) q[1];
sx q[1];
rz(-0.95744915) q[1];
rz(-2.0285151) q[3];
sx q[3];
rz(-1.7594171) q[3];
sx q[3];
rz(-0.45841226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6371969) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(1.8782328) q[2];
rz(1.0697621) q[3];
sx q[3];
rz(-1.0354778) q[3];
sx q[3];
rz(0.076816471) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99061154) q[0];
sx q[0];
rz(-0.69775109) q[0];
sx q[0];
rz(1.3569111) q[0];
rz(1.4487343) q[1];
sx q[1];
rz(-2.3285463) q[1];
sx q[1];
rz(2.9978233) q[1];
rz(2.8303296) q[2];
sx q[2];
rz(-2.0743721) q[2];
sx q[2];
rz(2.5596361) q[2];
rz(2.1830758) q[3];
sx q[3];
rz(-0.41194852) q[3];
sx q[3];
rz(1.3923399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
