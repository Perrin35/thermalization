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
rz(-1.5675867) q[1];
sx q[1];
rz(6.114638) q[1];
sx q[1];
rz(10.002879) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12500873) q[0];
sx q[0];
rz(-1.8468231) q[0];
sx q[0];
rz(-2.2651423) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69845818) q[2];
sx q[2];
rz(-1.466481) q[2];
sx q[2];
rz(-1.8836752) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4207124) q[1];
sx q[1];
rz(-1.136275) q[1];
sx q[1];
rz(-2.0220533) q[1];
x q[2];
rz(-0.16584478) q[3];
sx q[3];
rz(-0.87393099) q[3];
sx q[3];
rz(0.50717411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1987004) q[2];
sx q[2];
rz(-0.41112021) q[2];
sx q[2];
rz(1.6655507) q[2];
rz(0.57925159) q[3];
sx q[3];
rz(-1.9871291) q[3];
sx q[3];
rz(-0.66592413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85775527) q[0];
sx q[0];
rz(-1.0455766) q[0];
sx q[0];
rz(-0.055140821) q[0];
rz(2.227123) q[1];
sx q[1];
rz(-1.6998467) q[1];
sx q[1];
rz(1.2158016) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8912635) q[0];
sx q[0];
rz(-3.0760652) q[0];
sx q[0];
rz(0.28732462) q[0];
x q[1];
rz(-0.52771215) q[2];
sx q[2];
rz(-1.1368504) q[2];
sx q[2];
rz(2.3134856) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81406528) q[1];
sx q[1];
rz(-1.1720997) q[1];
sx q[1];
rz(-1.5908123) q[1];
rz(0.27894809) q[3];
sx q[3];
rz(-2.6761645) q[3];
sx q[3];
rz(2.5347749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2083644) q[2];
sx q[2];
rz(-2.6231982) q[2];
sx q[2];
rz(-2.5210157) q[2];
rz(-1.6879451) q[3];
sx q[3];
rz(-2.1098638) q[3];
sx q[3];
rz(-2.0645963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2700014) q[0];
sx q[0];
rz(-2.8102165) q[0];
sx q[0];
rz(1.6312067) q[0];
rz(0.67391467) q[1];
sx q[1];
rz(-1.8464512) q[1];
sx q[1];
rz(-1.5162226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1103863) q[0];
sx q[0];
rz(-1.0403311) q[0];
sx q[0];
rz(0.028859617) q[0];
rz(-pi) q[1];
rz(0.38336945) q[2];
sx q[2];
rz(-2.1738833) q[2];
sx q[2];
rz(2.4281339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5235345) q[1];
sx q[1];
rz(-2.6574572) q[1];
sx q[1];
rz(0.50393288) q[1];
rz(-0.45511873) q[3];
sx q[3];
rz(-0.60448217) q[3];
sx q[3];
rz(-2.4569619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7565833) q[2];
sx q[2];
rz(-1.5417121) q[2];
sx q[2];
rz(-0.6655244) q[2];
rz(-0.74337983) q[3];
sx q[3];
rz(-2.0205108) q[3];
sx q[3];
rz(0.22751787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3709563) q[0];
sx q[0];
rz(-1.0524858) q[0];
sx q[0];
rz(1.3002522) q[0];
rz(0.11996732) q[1];
sx q[1];
rz(-1.080546) q[1];
sx q[1];
rz(2.0948476) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20489731) q[0];
sx q[0];
rz(-3.129382) q[0];
sx q[0];
rz(-1.3365251) q[0];
x q[1];
rz(-1.98498) q[2];
sx q[2];
rz(-2.3799172) q[2];
sx q[2];
rz(-1.7047395) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0485301) q[1];
sx q[1];
rz(-1.0849864) q[1];
sx q[1];
rz(3.0792774) q[1];
x q[2];
rz(-2.1943521) q[3];
sx q[3];
rz(-0.6886607) q[3];
sx q[3];
rz(-0.50834828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3469424) q[2];
sx q[2];
rz(-2.0796937) q[2];
sx q[2];
rz(0.77155716) q[2];
rz(1.0509342) q[3];
sx q[3];
rz(-2.1080878) q[3];
sx q[3];
rz(1.1933901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0070888) q[0];
sx q[0];
rz(-1.6935885) q[0];
sx q[0];
rz(0.80528468) q[0];
rz(-1.6979506) q[1];
sx q[1];
rz(-2.3256358) q[1];
sx q[1];
rz(3.1051292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591025) q[0];
sx q[0];
rz(-0.56942372) q[0];
sx q[0];
rz(1.6605366) q[0];
x q[1];
rz(-1.6777322) q[2];
sx q[2];
rz(-0.57079878) q[2];
sx q[2];
rz(-2.609848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.441022) q[1];
sx q[1];
rz(-0.4800969) q[1];
sx q[1];
rz(-1.455485) q[1];
rz(-pi) q[2];
rz(1.7355326) q[3];
sx q[3];
rz(-2.142557) q[3];
sx q[3];
rz(-2.4248707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5263851) q[2];
sx q[2];
rz(-2.0168763) q[2];
sx q[2];
rz(0.57662326) q[2];
rz(1.5455101) q[3];
sx q[3];
rz(-2.2871064) q[3];
sx q[3];
rz(0.57687783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1103519) q[0];
sx q[0];
rz(-1.4875655) q[0];
sx q[0];
rz(-0.59979576) q[0];
rz(0.80232969) q[1];
sx q[1];
rz(-2.0498514) q[1];
sx q[1];
rz(-2.4937627) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9947056) q[0];
sx q[0];
rz(-2.0118666) q[0];
sx q[0];
rz(0.85710454) q[0];
x q[1];
rz(-3.134583) q[2];
sx q[2];
rz(-1.3607549) q[2];
sx q[2];
rz(1.5782859) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.036025612) q[1];
sx q[1];
rz(-0.67076761) q[1];
sx q[1];
rz(2.6545866) q[1];
x q[2];
rz(0.92354539) q[3];
sx q[3];
rz(-0.97935533) q[3];
sx q[3];
rz(-0.54852911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1442147) q[2];
sx q[2];
rz(-2.7684559) q[2];
sx q[2];
rz(-1.0972265) q[2];
rz(-2.1413474) q[3];
sx q[3];
rz(-2.2179243) q[3];
sx q[3];
rz(-1.3057115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10709396) q[0];
sx q[0];
rz(-1.9871563) q[0];
sx q[0];
rz(2.9423998) q[0];
rz(1.2983324) q[1];
sx q[1];
rz(-0.36111626) q[1];
sx q[1];
rz(0.64204204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605153) q[0];
sx q[0];
rz(-2.2937728) q[0];
sx q[0];
rz(-2.8882746) q[0];
rz(-pi) q[1];
rz(2.726905) q[2];
sx q[2];
rz(-1.0806335) q[2];
sx q[2];
rz(1.3140334) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.047028001) q[1];
sx q[1];
rz(-0.21512261) q[1];
sx q[1];
rz(2.7568629) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1690833) q[3];
sx q[3];
rz(-2.6582432) q[3];
sx q[3];
rz(-0.95248896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9127427) q[2];
sx q[2];
rz(-1.0656837) q[2];
sx q[2];
rz(1.9672811) q[2];
rz(-2.2187388) q[3];
sx q[3];
rz(-2.703981) q[3];
sx q[3];
rz(-2.9722884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4787503) q[0];
sx q[0];
rz(-1.4816544) q[0];
sx q[0];
rz(-0.99916512) q[0];
rz(-0.83813465) q[1];
sx q[1];
rz(-2.2331388) q[1];
sx q[1];
rz(2.1160486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89436707) q[0];
sx q[0];
rz(-2.0400088) q[0];
sx q[0];
rz(-0.85459015) q[0];
rz(-pi) q[1];
rz(-0.77825688) q[2];
sx q[2];
rz(-2.7822251) q[2];
sx q[2];
rz(3.0148413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9165695) q[1];
sx q[1];
rz(-0.65316641) q[1];
sx q[1];
rz(-2.2520425) q[1];
x q[2];
rz(2.4225967) q[3];
sx q[3];
rz(-1.1482802) q[3];
sx q[3];
rz(0.36976323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.020236882) q[2];
sx q[2];
rz(-1.3573703) q[2];
sx q[2];
rz(-0.21719246) q[2];
rz(-1.1413261) q[3];
sx q[3];
rz(-0.3749899) q[3];
sx q[3];
rz(-0.91514897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4311669) q[0];
sx q[0];
rz(-1.2942261) q[0];
sx q[0];
rz(-0.024209484) q[0];
rz(-1.0912033) q[1];
sx q[1];
rz(-1.1187436) q[1];
sx q[1];
rz(-2.3275163) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7221927) q[0];
sx q[0];
rz(-2.2792313) q[0];
sx q[0];
rz(-1.4577775) q[0];
rz(2.1948959) q[2];
sx q[2];
rz(-1.1128508) q[2];
sx q[2];
rz(0.68489546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7634552) q[1];
sx q[1];
rz(-2.3296851) q[1];
sx q[1];
rz(0.0017536963) q[1];
rz(0.086643593) q[3];
sx q[3];
rz(-1.3431755) q[3];
sx q[3];
rz(1.8342474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39390627) q[2];
sx q[2];
rz(-2.2795491) q[2];
sx q[2];
rz(1.5300592) q[2];
rz(-1.0511506) q[3];
sx q[3];
rz(-1.3408778) q[3];
sx q[3];
rz(0.10147258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53973389) q[0];
sx q[0];
rz(-2.1506385) q[0];
sx q[0];
rz(-0.47716004) q[0];
rz(1.5902144) q[1];
sx q[1];
rz(-1.662622) q[1];
sx q[1];
rz(-1.8345376) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73471776) q[0];
sx q[0];
rz(-0.85556036) q[0];
sx q[0];
rz(3.0815794) q[0];
rz(-1.1427504) q[2];
sx q[2];
rz(-1.1944886) q[2];
sx q[2];
rz(1.0264645) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.160678) q[1];
sx q[1];
rz(-2.0956934) q[1];
sx q[1];
rz(-0.69634931) q[1];
rz(-0.99193345) q[3];
sx q[3];
rz(-2.4483878) q[3];
sx q[3];
rz(-2.7424911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16605475) q[2];
sx q[2];
rz(-2.1259978) q[2];
sx q[2];
rz(-0.75237742) q[2];
rz(0.11263975) q[3];
sx q[3];
rz(-0.17387667) q[3];
sx q[3];
rz(-1.1627452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1174711) q[0];
sx q[0];
rz(-1.385067) q[0];
sx q[0];
rz(-2.0140482) q[0];
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
rz(2.9259089) q[3];
sx q[3];
rz(-0.75906039) q[3];
sx q[3];
rz(-1.0495521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
