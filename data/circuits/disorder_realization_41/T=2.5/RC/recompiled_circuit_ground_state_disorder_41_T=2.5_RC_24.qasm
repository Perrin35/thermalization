OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(1.289225) q[0];
sx q[0];
rz(9.3318648) q[0];
rz(3.0463123) q[1];
sx q[1];
rz(-2.4089101) q[1];
sx q[1];
rz(-1.9021775) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5372972) q[0];
sx q[0];
rz(-1.822246) q[0];
sx q[0];
rz(-2.9049113) q[0];
rz(2.8777166) q[2];
sx q[2];
rz(-1.3812764) q[2];
sx q[2];
rz(-0.60631982) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69529205) q[1];
sx q[1];
rz(-1.9009034) q[1];
sx q[1];
rz(-1.6020726) q[1];
rz(2.4437583) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(0.35240155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6713082) q[2];
sx q[2];
rz(-1.4802063) q[2];
sx q[2];
rz(-1.3489464) q[2];
rz(-2.3979483) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(2.9644137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0002366) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(0.29362383) q[0];
rz(-1.0307505) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(-2.7986599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1200877) q[0];
sx q[0];
rz(-1.7762799) q[0];
sx q[0];
rz(-2.1028825) q[0];
rz(-pi) q[1];
rz(1.2478825) q[2];
sx q[2];
rz(-0.39034778) q[2];
sx q[2];
rz(0.54803145) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82417497) q[1];
sx q[1];
rz(-1.7848099) q[1];
sx q[1];
rz(-0.70862464) q[1];
rz(0.94128709) q[3];
sx q[3];
rz(-1.3664748) q[3];
sx q[3];
rz(-0.91703892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12884101) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(-2.3823605) q[2];
rz(-0.020708474) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(-2.7032963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.52221209) q[0];
sx q[0];
rz(-0.601957) q[0];
sx q[0];
rz(0.0044599175) q[0];
rz(2.4941173) q[1];
sx q[1];
rz(-0.59440333) q[1];
sx q[1];
rz(2.3562145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7145536) q[0];
sx q[0];
rz(-1.6511962) q[0];
sx q[0];
rz(2.9820082) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0702281) q[2];
sx q[2];
rz(-2.7874261) q[2];
sx q[2];
rz(-3.0234697) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33357802) q[1];
sx q[1];
rz(-1.6393264) q[1];
sx q[1];
rz(1.0295111) q[1];
rz(2.249751) q[3];
sx q[3];
rz(-2.1168609) q[3];
sx q[3];
rz(-2.8308656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7711827) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.7837589) q[2];
rz(-0.34058288) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(-0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1444645) q[0];
sx q[0];
rz(-1.0839394) q[0];
sx q[0];
rz(0.88515627) q[0];
rz(-2.3947233) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(2.0451827) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0521419) q[0];
sx q[0];
rz(-2.6656287) q[0];
sx q[0];
rz(-0.83849283) q[0];
rz(-1.6302413) q[2];
sx q[2];
rz(-1.5402093) q[2];
sx q[2];
rz(-0.93418834) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8453522) q[1];
sx q[1];
rz(-1.8295604) q[1];
sx q[1];
rz(-2.6654763) q[1];
rz(-pi) q[2];
rz(-0.55620749) q[3];
sx q[3];
rz(-0.81327754) q[3];
sx q[3];
rz(2.5089056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(3.1169685) q[2];
rz(-2.5060182) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-0.88328254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6898952) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(0.46689335) q[0];
rz(3.0310071) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(-1.0708403) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2617874) q[0];
sx q[0];
rz(-0.51934411) q[0];
sx q[0];
rz(2.065361) q[0];
rz(-pi) q[1];
rz(-1.7532639) q[2];
sx q[2];
rz(-1.8844205) q[2];
sx q[2];
rz(-2.8387866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5683978) q[1];
sx q[1];
rz(-2.0802976) q[1];
sx q[1];
rz(2.0054818) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91584622) q[3];
sx q[3];
rz(-2.1612242) q[3];
sx q[3];
rz(2.6797323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0232627) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(0.037633745) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(-0.5948624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1325876) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(1.5124849) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(1.1423133) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61296755) q[0];
sx q[0];
rz(-2.7745947) q[0];
sx q[0];
rz(2.8662843) q[0];
rz(-pi) q[1];
rz(-3.0794607) q[2];
sx q[2];
rz(-1.3634063) q[2];
sx q[2];
rz(-1.0203938) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76872494) q[1];
sx q[1];
rz(-1.3228251) q[1];
sx q[1];
rz(-0.3526814) q[1];
rz(-pi) q[2];
rz(2.6352245) q[3];
sx q[3];
rz(-0.35249235) q[3];
sx q[3];
rz(1.6329721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3219354) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(-3.0211871) q[2];
rz(1.985792) q[3];
sx q[3];
rz(-1.5294231) q[3];
sx q[3];
rz(2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8026546) q[0];
sx q[0];
rz(-1.3405565) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(-1.5441719) q[1];
sx q[1];
rz(-1.495196) q[1];
sx q[1];
rz(-2.6905751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7190349) q[0];
sx q[0];
rz(-1.2698049) q[0];
sx q[0];
rz(1.6331571) q[0];
rz(-pi) q[1];
rz(-2.1233929) q[2];
sx q[2];
rz(-2.7180053) q[2];
sx q[2];
rz(-2.6747963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.800182) q[1];
sx q[1];
rz(-1.6918315) q[1];
sx q[1];
rz(0.43612632) q[1];
x q[2];
rz(-2.3312206) q[3];
sx q[3];
rz(-1.5640508) q[3];
sx q[3];
rz(-0.41616671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94577998) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(-1.7808524) q[2];
rz(-2.1360548) q[3];
sx q[3];
rz(-1.1755627) q[3];
sx q[3];
rz(-1.7520693) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239531) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(0.70607591) q[0];
rz(1.5977244) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(0.85618883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2285952) q[0];
sx q[0];
rz(-1.5902336) q[0];
sx q[0];
rz(0.0073846505) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43120877) q[2];
sx q[2];
rz(-2.2715306) q[2];
sx q[2];
rz(-1.7675811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.013102839) q[1];
sx q[1];
rz(-1.0016514) q[1];
sx q[1];
rz(-1.8123958) q[1];
rz(-pi) q[2];
rz(1.0821277) q[3];
sx q[3];
rz(-1.6889945) q[3];
sx q[3];
rz(1.7080054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2989444) q[2];
sx q[2];
rz(-1.7405258) q[2];
sx q[2];
rz(0.16743463) q[2];
rz(-1.5873448) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(1.0003482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.7635968) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(-0.3717306) q[0];
rz(1.7077712) q[1];
sx q[1];
rz(-1.6038409) q[1];
sx q[1];
rz(2.8415714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2903295) q[0];
sx q[0];
rz(-1.4545049) q[0];
sx q[0];
rz(-0.28684464) q[0];
rz(-pi) q[1];
rz(-1.5302363) q[2];
sx q[2];
rz(-2.1475361) q[2];
sx q[2];
rz(-0.64740136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65574291) q[1];
sx q[1];
rz(-0.92461899) q[1];
sx q[1];
rz(1.7263401) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12436538) q[3];
sx q[3];
rz(-0.87166407) q[3];
sx q[3];
rz(-1.3209526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6649449) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(0.42759582) q[2];
rz(1.556373) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(0.75470406) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4204191) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(-2.6721201) q[0];
rz(0.92533127) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(3.1184149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0231004) q[0];
sx q[0];
rz(-1.4064624) q[0];
sx q[0];
rz(-2.611155) q[0];
rz(2.2830354) q[2];
sx q[2];
rz(-2.3344667) q[2];
sx q[2];
rz(-1.6940534) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1851481) q[1];
sx q[1];
rz(-1.4063764) q[1];
sx q[1];
rz(-1.7471353) q[1];
rz(1.7816824) q[3];
sx q[3];
rz(-0.74112219) q[3];
sx q[3];
rz(-1.1531342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8699441) q[2];
sx q[2];
rz(-1.9359438) q[2];
sx q[2];
rz(-0.49087697) q[2];
rz(1.0902181) q[3];
sx q[3];
rz(-0.96929437) q[3];
sx q[3];
rz(0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.28521095) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(-0.27676997) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(2.7748952) q[2];
sx q[2];
rz(-0.72281217) q[2];
sx q[2];
rz(-2.6642961) q[2];
rz(-3.1076486) q[3];
sx q[3];
rz(-2.6068519) q[3];
sx q[3];
rz(3.1153771) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
