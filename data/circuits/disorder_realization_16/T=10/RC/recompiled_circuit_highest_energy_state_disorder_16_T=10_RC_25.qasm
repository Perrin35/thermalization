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
rz(-1.1666522) q[0];
sx q[0];
rz(-0.44297543) q[0];
sx q[0];
rz(-0.38323453) q[0];
rz(-1.2869599) q[1];
sx q[1];
rz(-1.1416963) q[1];
sx q[1];
rz(-0.72007522) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0135277) q[0];
sx q[0];
rz(-1.9682471) q[0];
sx q[0];
rz(1.1532408) q[0];
rz(-pi) q[1];
x q[1];
rz(1.371648) q[2];
sx q[2];
rz(-1.3825644) q[2];
sx q[2];
rz(-0.07478274) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.076700828) q[1];
sx q[1];
rz(-1.2982839) q[1];
sx q[1];
rz(0.38857675) q[1];
rz(0.96020697) q[3];
sx q[3];
rz(-1.4355324) q[3];
sx q[3];
rz(-1.3162515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2893452) q[2];
sx q[2];
rz(-0.87163681) q[2];
sx q[2];
rz(2.3343425) q[2];
rz(1.4489669) q[3];
sx q[3];
rz(-0.92985409) q[3];
sx q[3];
rz(-0.71678954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.1194696) q[0];
sx q[0];
rz(-1.2998281) q[0];
sx q[0];
rz(-1.8362554) q[0];
rz(-2.3985825) q[1];
sx q[1];
rz(-0.65736714) q[1];
sx q[1];
rz(-1.8441127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3548879) q[0];
sx q[0];
rz(-0.44092766) q[0];
sx q[0];
rz(-1.132325) q[0];
rz(-pi) q[1];
rz(2.2114685) q[2];
sx q[2];
rz(-0.67652297) q[2];
sx q[2];
rz(-0.97201024) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.686649) q[1];
sx q[1];
rz(-1.8669858) q[1];
sx q[1];
rz(0.95570968) q[1];
rz(2.7801974) q[3];
sx q[3];
rz(-1.6797429) q[3];
sx q[3];
rz(-1.7360842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.3798736) q[2];
sx q[2];
rz(-0.97731176) q[2];
sx q[2];
rz(1.0404111) q[2];
rz(-2.8271293) q[3];
sx q[3];
rz(-1.9324666) q[3];
sx q[3];
rz(1.2578806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8085025) q[0];
sx q[0];
rz(-2.1603778) q[0];
sx q[0];
rz(1.8044949) q[0];
rz(2.5143276) q[1];
sx q[1];
rz(-1.4074872) q[1];
sx q[1];
rz(-1.6076535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8117789) q[0];
sx q[0];
rz(-2.3372991) q[0];
sx q[0];
rz(0.95383184) q[0];
rz(-pi) q[1];
rz(-1.8684325) q[2];
sx q[2];
rz(-1.9797851) q[2];
sx q[2];
rz(1.2163509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7067553) q[1];
sx q[1];
rz(-1.1011775) q[1];
sx q[1];
rz(-1.7781185) q[1];
rz(-2.0407143) q[3];
sx q[3];
rz(-0.4199314) q[3];
sx q[3];
rz(0.41929276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4072121) q[2];
sx q[2];
rz(-2.247159) q[2];
sx q[2];
rz(1.2073995) q[2];
rz(-2.1842128) q[3];
sx q[3];
rz(-1.9060241) q[3];
sx q[3];
rz(-0.68937075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44416881) q[0];
sx q[0];
rz(-2.1902695) q[0];
sx q[0];
rz(0.091766894) q[0];
rz(0.16608876) q[1];
sx q[1];
rz(-1.2867915) q[1];
sx q[1];
rz(2.1002358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9516141) q[0];
sx q[0];
rz(-2.3654177) q[0];
sx q[0];
rz(-1.6329732) q[0];
rz(-1.1686334) q[2];
sx q[2];
rz(-0.82457525) q[2];
sx q[2];
rz(-2.4006725) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68616679) q[1];
sx q[1];
rz(-2.0435026) q[1];
sx q[1];
rz(3.1030639) q[1];
rz(2.5934173) q[3];
sx q[3];
rz(-2.976307) q[3];
sx q[3];
rz(2.838244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.912821) q[2];
sx q[2];
rz(-0.81835881) q[2];
sx q[2];
rz(-1.5405601) q[2];
rz(2.1238964) q[3];
sx q[3];
rz(-1.8120268) q[3];
sx q[3];
rz(-1.8890007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8215264) q[0];
sx q[0];
rz(-1.3952661) q[0];
sx q[0];
rz(-0.23859247) q[0];
rz(2.353031) q[1];
sx q[1];
rz(-2.8193654) q[1];
sx q[1];
rz(2.2408748) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.032388587) q[0];
sx q[0];
rz(-0.51133832) q[0];
sx q[0];
rz(2.5912259) q[0];
rz(-1.3447755) q[2];
sx q[2];
rz(-0.39835762) q[2];
sx q[2];
rz(-0.67581723) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5594192) q[1];
sx q[1];
rz(-1.6985849) q[1];
sx q[1];
rz(-0.41499563) q[1];
rz(-2.4191148) q[3];
sx q[3];
rz(-0.9831419) q[3];
sx q[3];
rz(2.8624299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0476524) q[2];
sx q[2];
rz(-1.36146) q[2];
sx q[2];
rz(1.368604) q[2];
rz(0.37745825) q[3];
sx q[3];
rz(-2.3182175) q[3];
sx q[3];
rz(-1.1641758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222328) q[0];
sx q[0];
rz(-0.43411175) q[0];
sx q[0];
rz(2.9121616) q[0];
rz(-2.5550628) q[1];
sx q[1];
rz(-2.6507381) q[1];
sx q[1];
rz(1.4803001) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.081588) q[0];
sx q[0];
rz(-1.2179273) q[0];
sx q[0];
rz(0.88624432) q[0];
rz(-pi) q[1];
rz(2.7469559) q[2];
sx q[2];
rz(-0.60606128) q[2];
sx q[2];
rz(0.74575033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6203998) q[1];
sx q[1];
rz(-1.344465) q[1];
sx q[1];
rz(1.625555) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.839572) q[3];
sx q[3];
rz(-2.44938) q[3];
sx q[3];
rz(0.48066329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71969879) q[2];
sx q[2];
rz(-0.59013683) q[2];
sx q[2];
rz(-1.3479007) q[2];
rz(0.19139309) q[3];
sx q[3];
rz(-2.8189711) q[3];
sx q[3];
rz(1.9026683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.557068) q[0];
sx q[0];
rz(-1.6275591) q[0];
sx q[0];
rz(0.53644449) q[0];
rz(0.084511936) q[1];
sx q[1];
rz(-2.5630496) q[1];
sx q[1];
rz(1.1892148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3670334) q[0];
sx q[0];
rz(-1.7248472) q[0];
sx q[0];
rz(0.94512248) q[0];
x q[1];
rz(0.68185735) q[2];
sx q[2];
rz(-1.8769662) q[2];
sx q[2];
rz(0.078753565) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8457905) q[1];
sx q[1];
rz(-1.2582878) q[1];
sx q[1];
rz(2.0927696) q[1];
rz(-pi) q[2];
rz(1.1105082) q[3];
sx q[3];
rz(-1.7841859) q[3];
sx q[3];
rz(1.8387356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9540003) q[2];
sx q[2];
rz(-1.2728929) q[2];
sx q[2];
rz(-3.0583337) q[2];
rz(1.8795053) q[3];
sx q[3];
rz(-0.83674651) q[3];
sx q[3];
rz(-0.73778233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5872203) q[0];
sx q[0];
rz(-1.6393336) q[0];
sx q[0];
rz(-0.7837522) q[0];
rz(1.6571677) q[1];
sx q[1];
rz(-0.52394167) q[1];
sx q[1];
rz(-0.2934244) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992511) q[0];
sx q[0];
rz(-0.60428491) q[0];
sx q[0];
rz(2.2544331) q[0];
rz(-0.22343724) q[2];
sx q[2];
rz(-0.91470892) q[2];
sx q[2];
rz(-0.81435299) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3345016) q[1];
sx q[1];
rz(-1.1576022) q[1];
sx q[1];
rz(-2.9838802) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.226172) q[3];
sx q[3];
rz(-1.3815945) q[3];
sx q[3];
rz(-0.20841852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36711127) q[2];
sx q[2];
rz(-1.7256871) q[2];
sx q[2];
rz(-1.6488546) q[2];
rz(-1.1501009) q[3];
sx q[3];
rz(-0.29905683) q[3];
sx q[3];
rz(-2.2404631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8915326) q[0];
sx q[0];
rz(-2.1136668) q[0];
sx q[0];
rz(0.0082536396) q[0];
rz(0.0025657733) q[1];
sx q[1];
rz(-0.51104128) q[1];
sx q[1];
rz(-1.0645688) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24495793) q[0];
sx q[0];
rz(-1.1616529) q[0];
sx q[0];
rz(-0.70434241) q[0];
rz(1.538689) q[2];
sx q[2];
rz(-0.24954441) q[2];
sx q[2];
rz(-1.3448993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4618116) q[1];
sx q[1];
rz(-1.0625328) q[1];
sx q[1];
rz(-2.630788) q[1];
rz(-2.921815) q[3];
sx q[3];
rz(-0.6005601) q[3];
sx q[3];
rz(-1.8547577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.061444) q[2];
sx q[2];
rz(-1.6950357) q[2];
sx q[2];
rz(-0.20305571) q[2];
rz(-1.9301346) q[3];
sx q[3];
rz(-1.0695499) q[3];
sx q[3];
rz(-3.0058461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1900629) q[0];
sx q[0];
rz(-2.0385346) q[0];
sx q[0];
rz(2.4679825) q[0];
rz(-2.8296962) q[1];
sx q[1];
rz(-1.936828) q[1];
sx q[1];
rz(-1.9685251) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86248428) q[0];
sx q[0];
rz(-1.0432093) q[0];
sx q[0];
rz(-0.34054813) q[0];
rz(-2.8946213) q[2];
sx q[2];
rz(-0.77148831) q[2];
sx q[2];
rz(2.889973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0153326) q[1];
sx q[1];
rz(-1.4950206) q[1];
sx q[1];
rz(-1.9858236) q[1];
rz(-2.2820149) q[3];
sx q[3];
rz(-1.2378581) q[3];
sx q[3];
rz(-0.22151079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6761026) q[2];
sx q[2];
rz(-0.61773053) q[2];
sx q[2];
rz(-0.22107302) q[2];
rz(-0.57340932) q[3];
sx q[3];
rz(-0.87380496) q[3];
sx q[3];
rz(-1.1055841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318646) q[0];
sx q[0];
rz(-1.8479713) q[0];
sx q[0];
rz(0.027298409) q[0];
rz(-1.191054) q[1];
sx q[1];
rz(-1.4203032) q[1];
sx q[1];
rz(1.3994093) q[1];
rz(-1.9223735) q[2];
sx q[2];
rz(-0.58813358) q[2];
sx q[2];
rz(1.0434601) q[2];
rz(2.7577362) q[3];
sx q[3];
rz(-1.0105269) q[3];
sx q[3];
rz(-0.15906048) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
