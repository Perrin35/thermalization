OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2616413) q[0];
sx q[0];
rz(-0.14544848) q[0];
sx q[0];
rz(-2.8737336) q[0];
rz(1.574006) q[1];
sx q[1];
rz(-2.9730453) q[1];
sx q[1];
rz(-0.57810098) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12500873) q[0];
sx q[0];
rz(-1.2947695) q[0];
sx q[0];
rz(0.87645032) q[0];
rz(-pi) q[1];
rz(0.16139754) q[2];
sx q[2];
rz(-0.70490743) q[2];
sx q[2];
rz(-2.7053331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0069258) q[1];
sx q[1];
rz(-2.5257266) q[1];
sx q[1];
rz(-2.3872972) q[1];
x q[2];
rz(-1.7655444) q[3];
sx q[3];
rz(-2.4284902) q[3];
sx q[3];
rz(-2.3793067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.1987004) q[2];
sx q[2];
rz(-2.7304724) q[2];
sx q[2];
rz(1.6655507) q[2];
rz(0.57925159) q[3];
sx q[3];
rz(-1.1544635) q[3];
sx q[3];
rz(-2.4756685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85775527) q[0];
sx q[0];
rz(-1.0455766) q[0];
sx q[0];
rz(0.055140821) q[0];
rz(2.227123) q[1];
sx q[1];
rz(-1.6998467) q[1];
sx q[1];
rz(1.2158016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8912635) q[0];
sx q[0];
rz(-3.0760652) q[0];
sx q[0];
rz(2.854268) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74389768) q[2];
sx q[2];
rz(-2.4716931) q[2];
sx q[2];
rz(0.11775859) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86558531) q[1];
sx q[1];
rz(-2.7424208) q[1];
sx q[1];
rz(-0.047476032) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8626446) q[3];
sx q[3];
rz(-2.6761645) q[3];
sx q[3];
rz(0.60681776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9332283) q[2];
sx q[2];
rz(-0.51839447) q[2];
sx q[2];
rz(-0.62057692) q[2];
rz(-1.4536475) q[3];
sx q[3];
rz(-2.1098638) q[3];
sx q[3];
rz(2.0645963) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2700014) q[0];
sx q[0];
rz(-2.8102165) q[0];
sx q[0];
rz(-1.6312067) q[0];
rz(-2.467678) q[1];
sx q[1];
rz(-1.2951415) q[1];
sx q[1];
rz(-1.6253701) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1103863) q[0];
sx q[0];
rz(-2.1012615) q[0];
sx q[0];
rz(-0.028859617) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2095334) q[2];
sx q[2];
rz(-1.2576418) q[2];
sx q[2];
rz(-2.0593638) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4067799) q[1];
sx q[1];
rz(-1.3441097) q[1];
sx q[1];
rz(-0.43155687) q[1];
x q[2];
rz(1.8655769) q[3];
sx q[3];
rz(-1.0350482) q[3];
sx q[3];
rz(1.2212041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7565833) q[2];
sx q[2];
rz(-1.5417121) q[2];
sx q[2];
rz(0.6655244) q[2];
rz(-0.74337983) q[3];
sx q[3];
rz(-1.1210818) q[3];
sx q[3];
rz(-0.22751787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7706364) q[0];
sx q[0];
rz(-1.0524858) q[0];
sx q[0];
rz(1.8413405) q[0];
rz(0.11996732) q[1];
sx q[1];
rz(-2.0610466) q[1];
sx q[1];
rz(1.0467451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20489731) q[0];
sx q[0];
rz(-3.129382) q[0];
sx q[0];
rz(-1.3365251) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1566126) q[2];
sx q[2];
rz(-2.3799172) q[2];
sx q[2];
rz(1.7047395) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6347329) q[1];
sx q[1];
rz(-1.5156989) q[1];
sx q[1];
rz(2.0574089) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1943521) q[3];
sx q[3];
rz(-0.6886607) q[3];
sx q[3];
rz(0.50834828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3469424) q[2];
sx q[2];
rz(-2.0796937) q[2];
sx q[2];
rz(-2.3700355) q[2];
rz(2.0906585) q[3];
sx q[3];
rz(-1.0335048) q[3];
sx q[3];
rz(1.1933901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0070888) q[0];
sx q[0];
rz(-1.4480042) q[0];
sx q[0];
rz(-0.80528468) q[0];
rz(1.6979506) q[1];
sx q[1];
rz(-0.81595683) q[1];
sx q[1];
rz(3.1051292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5492461) q[0];
sx q[0];
rz(-2.1376462) q[0];
sx q[0];
rz(3.0842848) q[0];
rz(-pi) q[1];
rz(-1.4638605) q[2];
sx q[2];
rz(-0.57079878) q[2];
sx q[2];
rz(2.609848) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.441022) q[1];
sx q[1];
rz(-0.4800969) q[1];
sx q[1];
rz(1.6861077) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5636118) q[3];
sx q[3];
rz(-1.4324463) q[3];
sx q[3];
rz(2.377233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61520758) q[2];
sx q[2];
rz(-2.0168763) q[2];
sx q[2];
rz(0.57662326) q[2];
rz(-1.5960826) q[3];
sx q[3];
rz(-2.2871064) q[3];
sx q[3];
rz(0.57687783) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0312408) q[0];
sx q[0];
rz(-1.6540271) q[0];
sx q[0];
rz(0.59979576) q[0];
rz(-2.339263) q[1];
sx q[1];
rz(-2.0498514) q[1];
sx q[1];
rz(0.64782992) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8819274) q[0];
sx q[0];
rz(-2.3234832) q[0];
sx q[0];
rz(2.1955793) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7808427) q[2];
sx q[2];
rz(-1.577652) q[2];
sx q[2];
rz(0.0060280212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.036025612) q[1];
sx q[1];
rz(-0.67076761) q[1];
sx q[1];
rz(-0.48700602) q[1];
x q[2];
rz(2.2180473) q[3];
sx q[3];
rz(-2.1622373) q[3];
sx q[3];
rz(2.5930635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1442147) q[2];
sx q[2];
rz(-0.37313676) q[2];
sx q[2];
rz(-2.0443661) q[2];
rz(-2.1413474) q[3];
sx q[3];
rz(-0.92366832) q[3];
sx q[3];
rz(1.3057115) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10709396) q[0];
sx q[0];
rz(-1.9871563) q[0];
sx q[0];
rz(-2.9423998) q[0];
rz(1.2983324) q[1];
sx q[1];
rz(-0.36111626) q[1];
sx q[1];
rz(0.64204204) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4810774) q[0];
sx q[0];
rz(-0.84781983) q[0];
sx q[0];
rz(-2.8882746) q[0];
rz(0.41468765) q[2];
sx q[2];
rz(-1.0806335) q[2];
sx q[2];
rz(-1.3140334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.047028001) q[1];
sx q[1];
rz(-2.92647) q[1];
sx q[1];
rz(2.7568629) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1690833) q[3];
sx q[3];
rz(-0.48334941) q[3];
sx q[3];
rz(2.1891037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9127427) q[2];
sx q[2];
rz(-1.0656837) q[2];
sx q[2];
rz(1.1743116) q[2];
rz(-2.2187388) q[3];
sx q[3];
rz(-0.43761161) q[3];
sx q[3];
rz(-0.16930425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4787503) q[0];
sx q[0];
rz(-1.4816544) q[0];
sx q[0];
rz(2.1424275) q[0];
rz(0.83813465) q[1];
sx q[1];
rz(-0.90845388) q[1];
sx q[1];
rz(-1.025544) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9860172) q[0];
sx q[0];
rz(-2.3087569) q[0];
sx q[0];
rz(2.2283594) q[0];
rz(-0.26141459) q[2];
sx q[2];
rz(-1.8202708) q[2];
sx q[2];
rz(0.69875137) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5708988) q[1];
sx q[1];
rz(-2.0624196) q[1];
sx q[1];
rz(2.6925283) q[1];
rz(-0.59896627) q[3];
sx q[3];
rz(-2.3271797) q[3];
sx q[3];
rz(1.6391476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.020236882) q[2];
sx q[2];
rz(-1.3573703) q[2];
sx q[2];
rz(2.9244002) q[2];
rz(-2.0002666) q[3];
sx q[3];
rz(-0.3749899) q[3];
sx q[3];
rz(-2.2264437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4311669) q[0];
sx q[0];
rz(-1.8473666) q[0];
sx q[0];
rz(3.1173832) q[0];
rz(-1.0912033) q[1];
sx q[1];
rz(-1.1187436) q[1];
sx q[1];
rz(-2.3275163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0776802) q[0];
sx q[0];
rz(-1.6565431) q[0];
sx q[0];
rz(2.4299939) q[0];
x q[1];
rz(-2.5957554) q[2];
sx q[2];
rz(-2.1225192) q[2];
sx q[2];
rz(-1.1941225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.947727) q[1];
sx q[1];
rz(-1.5695238) q[1];
sx q[1];
rz(-2.3296859) q[1];
rz(1.9283251) q[3];
sx q[3];
rz(-2.8983064) q[3];
sx q[3];
rz(1.4668087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7476864) q[2];
sx q[2];
rz(-2.2795491) q[2];
sx q[2];
rz(-1.6115335) q[2];
rz(2.090442) q[3];
sx q[3];
rz(-1.3408778) q[3];
sx q[3];
rz(0.10147258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6018588) q[0];
sx q[0];
rz(-0.99095416) q[0];
sx q[0];
rz(-0.47716004) q[0];
rz(-1.5902144) q[1];
sx q[1];
rz(-1.4789707) q[1];
sx q[1];
rz(-1.8345376) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982419) q[0];
sx q[0];
rz(-2.4242867) q[0];
sx q[0];
rz(1.5018612) q[0];
x q[1];
rz(-2.731852) q[2];
sx q[2];
rz(-1.9671408) q[2];
sx q[2];
rz(-0.71045638) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2714398) q[1];
sx q[1];
rz(-0.84476568) q[1];
sx q[1];
rz(2.4072263) q[1];
rz(-pi) q[2];
rz(-0.96316506) q[3];
sx q[3];
rz(-1.9279216) q[3];
sx q[3];
rz(-0.70588338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16605475) q[2];
sx q[2];
rz(-2.1259978) q[2];
sx q[2];
rz(-0.75237742) q[2];
rz(0.11263975) q[3];
sx q[3];
rz(-0.17387667) q[3];
sx q[3];
rz(1.9788474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.024121506) q[0];
sx q[0];
rz(-1.7565256) q[0];
sx q[0];
rz(1.1275445) q[0];
rz(-2.7680001) q[1];
sx q[1];
rz(-1.6289381) q[1];
sx q[1];
rz(-2.1388114) q[1];
rz(1.9890979) q[2];
sx q[2];
rz(-1.9807182) q[2];
sx q[2];
rz(2.9302927) q[2];
rz(1.3704902) q[3];
sx q[3];
rz(-2.3080993) q[3];
sx q[3];
rz(-0.75626683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
