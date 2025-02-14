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
rz(0.54230827) q[0];
sx q[0];
rz(-0.13442726) q[0];
sx q[0];
rz(-1.0472714) q[0];
rz(-0.32416999) q[1];
sx q[1];
rz(2.9919762) q[1];
sx q[1];
rz(11.621578) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1143998) q[0];
sx q[0];
rz(-1.8437705) q[0];
sx q[0];
rz(1.2370212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23350291) q[2];
sx q[2];
rz(-2.2036607) q[2];
sx q[2];
rz(-1.3695182) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45356645) q[1];
sx q[1];
rz(-2.0658148) q[1];
sx q[1];
rz(-1.051245) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7994491) q[3];
sx q[3];
rz(-0.49352577) q[3];
sx q[3];
rz(2.760104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85113877) q[2];
sx q[2];
rz(-1.6893427) q[2];
sx q[2];
rz(1.5645082) q[2];
rz(2.5751298) q[3];
sx q[3];
rz(-2.5201859) q[3];
sx q[3];
rz(-1.8433146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8241149) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(0.7199921) q[0];
rz(0.27436817) q[1];
sx q[1];
rz(-0.73148483) q[1];
sx q[1];
rz(2.5879587) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33687978) q[0];
sx q[0];
rz(-1.9229445) q[0];
sx q[0];
rz(3.0180305) q[0];
rz(-pi) q[1];
rz(0.64525147) q[2];
sx q[2];
rz(-0.49075365) q[2];
sx q[2];
rz(0.09749271) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8588176) q[1];
sx q[1];
rz(-1.2056279) q[1];
sx q[1];
rz(-2.3479985) q[1];
x q[2];
rz(-2.3838061) q[3];
sx q[3];
rz(-0.67000853) q[3];
sx q[3];
rz(-0.78950933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1729892) q[2];
sx q[2];
rz(-1.9393238) q[2];
sx q[2];
rz(0.36390641) q[2];
rz(0.11161741) q[3];
sx q[3];
rz(-0.93897396) q[3];
sx q[3];
rz(0.76044559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97742057) q[0];
sx q[0];
rz(-2.3168679) q[0];
sx q[0];
rz(0.0053996276) q[0];
rz(-0.81575704) q[1];
sx q[1];
rz(-1.1382256) q[1];
sx q[1];
rz(-0.38591787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621423) q[0];
sx q[0];
rz(-1.291692) q[0];
sx q[0];
rz(-0.35274679) q[0];
rz(2.1109525) q[2];
sx q[2];
rz(-0.6331501) q[2];
sx q[2];
rz(-1.6748705) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5623916) q[1];
sx q[1];
rz(-1.950197) q[1];
sx q[1];
rz(-3.0048278) q[1];
rz(-pi) q[2];
rz(1.8688275) q[3];
sx q[3];
rz(-2.7103817) q[3];
sx q[3];
rz(1.9995156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5526814) q[2];
sx q[2];
rz(-2.3083355) q[2];
sx q[2];
rz(1.7595278) q[2];
rz(-2.9803993) q[3];
sx q[3];
rz(-1.7101945) q[3];
sx q[3];
rz(2.5126357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88518077) q[0];
sx q[0];
rz(-1.8653402) q[0];
sx q[0];
rz(1.3389583) q[0];
rz(-1.8244052) q[1];
sx q[1];
rz(-2.1664797) q[1];
sx q[1];
rz(0.29979527) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6000344) q[0];
sx q[0];
rz(-2.2081714) q[0];
sx q[0];
rz(-1.1336898) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1607397) q[2];
sx q[2];
rz(-1.9741657) q[2];
sx q[2];
rz(1.174675) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2721709) q[1];
sx q[1];
rz(-1.4336839) q[1];
sx q[1];
rz(1.2161944) q[1];
rz(3.0092029) q[3];
sx q[3];
rz(-2.2023938) q[3];
sx q[3];
rz(-1.5204483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6224711) q[2];
sx q[2];
rz(-0.93834472) q[2];
sx q[2];
rz(0.31912121) q[2];
rz(-3.0506813) q[3];
sx q[3];
rz(-2.9967283) q[3];
sx q[3];
rz(1.3566141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262481) q[0];
sx q[0];
rz(-2.321796) q[0];
sx q[0];
rz(0.0023284624) q[0];
rz(-1.5709411) q[1];
sx q[1];
rz(-2.3410773) q[1];
sx q[1];
rz(1.9470107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1462458) q[0];
sx q[0];
rz(-1.6888535) q[0];
sx q[0];
rz(0.15583584) q[0];
x q[1];
rz(1.5244687) q[2];
sx q[2];
rz(-0.71832685) q[2];
sx q[2];
rz(2.3026932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4264226) q[1];
sx q[1];
rz(-1.914901) q[1];
sx q[1];
rz(-0.10870966) q[1];
rz(-pi) q[2];
rz(0.55641716) q[3];
sx q[3];
rz(-1.3297982) q[3];
sx q[3];
rz(-0.3410546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.883256) q[2];
sx q[2];
rz(-0.77815762) q[2];
sx q[2];
rz(0.077433132) q[2];
rz(-2.1954913) q[3];
sx q[3];
rz(-0.48282048) q[3];
sx q[3];
rz(-0.4536804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2133863) q[0];
sx q[0];
rz(-3.0551857) q[0];
sx q[0];
rz(1.7267831) q[0];
rz(0.66697031) q[1];
sx q[1];
rz(-1.357115) q[1];
sx q[1];
rz(-0.40474969) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0721829) q[0];
sx q[0];
rz(-1.2889922) q[0];
sx q[0];
rz(-2.826716) q[0];
rz(-pi) q[1];
rz(1.9046632) q[2];
sx q[2];
rz(-0.68511663) q[2];
sx q[2];
rz(1.8917055) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7915899) q[1];
sx q[1];
rz(-1.631893) q[1];
sx q[1];
rz(2.0339802) q[1];
x q[2];
rz(-2.0575298) q[3];
sx q[3];
rz(-1.2477418) q[3];
sx q[3];
rz(2.9655181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82379597) q[2];
sx q[2];
rz(-0.78936374) q[2];
sx q[2];
rz(0.44343597) q[2];
rz(-2.7277842) q[3];
sx q[3];
rz(-0.2897073) q[3];
sx q[3];
rz(2.2558291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39171788) q[0];
sx q[0];
rz(-1.3626008) q[0];
sx q[0];
rz(2.474127) q[0];
rz(-3.093847) q[1];
sx q[1];
rz(-1.1404488) q[1];
sx q[1];
rz(2.023229) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3781609) q[0];
sx q[0];
rz(-1.973828) q[0];
sx q[0];
rz(0.58149882) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1492858) q[2];
sx q[2];
rz(-1.4690781) q[2];
sx q[2];
rz(0.5660688) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0065919) q[1];
sx q[1];
rz(-1.9891095) q[1];
sx q[1];
rz(-2.5043284) q[1];
rz(-2.0703378) q[3];
sx q[3];
rz(-1.6304255) q[3];
sx q[3];
rz(-1.7891375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.095801) q[2];
sx q[2];
rz(-1.8818776) q[2];
sx q[2];
rz(-0.0087139159) q[2];
rz(1.9376612) q[3];
sx q[3];
rz(-0.34398505) q[3];
sx q[3];
rz(3.1194527) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1190204) q[0];
sx q[0];
rz(-2.75596) q[0];
sx q[0];
rz(2.2517396) q[0];
rz(-1.285137) q[1];
sx q[1];
rz(-1.0533918) q[1];
sx q[1];
rz(-0.13592517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3691007) q[0];
sx q[0];
rz(-1.6181503) q[0];
sx q[0];
rz(2.9124509) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68259062) q[2];
sx q[2];
rz(-2.3851378) q[2];
sx q[2];
rz(2.3025049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0057266) q[1];
sx q[1];
rz(-2.7512065) q[1];
sx q[1];
rz(-0.54061215) q[1];
x q[2];
rz(-2.9985524) q[3];
sx q[3];
rz(-1.2238127) q[3];
sx q[3];
rz(-0.59511371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4385684) q[2];
sx q[2];
rz(-2.2615304) q[2];
sx q[2];
rz(2.0795889) q[2];
rz(2.2117173) q[3];
sx q[3];
rz(-2.3558741) q[3];
sx q[3];
rz(2.3532531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7787665) q[0];
sx q[0];
rz(-2.7500948) q[0];
sx q[0];
rz(0.40089259) q[0];
rz(-0.76155424) q[1];
sx q[1];
rz(-1.6490033) q[1];
sx q[1];
rz(-0.78899312) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83037277) q[0];
sx q[0];
rz(-1.5867751) q[0];
sx q[0];
rz(-3.014084) q[0];
x q[1];
rz(0.15536552) q[2];
sx q[2];
rz(-0.90165686) q[2];
sx q[2];
rz(2.5071627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4915062) q[1];
sx q[1];
rz(-1.6900926) q[1];
sx q[1];
rz(-2.5507798) q[1];
rz(-pi) q[2];
rz(-0.085461334) q[3];
sx q[3];
rz(-0.52640593) q[3];
sx q[3];
rz(1.5861024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.57364982) q[2];
sx q[2];
rz(-1.8486134) q[2];
sx q[2];
rz(2.4207777) q[2];
rz(1.4912262) q[3];
sx q[3];
rz(-0.78250116) q[3];
sx q[3];
rz(-1.3097552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023130527) q[0];
sx q[0];
rz(-3.0266302) q[0];
sx q[0];
rz(-0.7777099) q[0];
rz(-0.65981162) q[1];
sx q[1];
rz(-0.86645627) q[1];
sx q[1];
rz(-2.7515817) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5725419) q[0];
sx q[0];
rz(-0.13992289) q[0];
sx q[0];
rz(1.8403017) q[0];
rz(-0.71163746) q[2];
sx q[2];
rz(-0.88274985) q[2];
sx q[2];
rz(1.9311063) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6811878) q[1];
sx q[1];
rz(-1.273566) q[1];
sx q[1];
rz(-1.1497208) q[1];
x q[2];
rz(0.94160752) q[3];
sx q[3];
rz(-2.3857085) q[3];
sx q[3];
rz(-0.98996491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4593792) q[2];
sx q[2];
rz(-2.7028658) q[2];
sx q[2];
rz(2.6518346) q[2];
rz(-0.30719906) q[3];
sx q[3];
rz(-0.92370737) q[3];
sx q[3];
rz(-1.3908305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3601892) q[0];
sx q[0];
rz(-1.6506945) q[0];
sx q[0];
rz(-1.6682464) q[0];
rz(-2.0669943) q[1];
sx q[1];
rz(-2.2886724) q[1];
sx q[1];
rz(1.2235175) q[1];
rz(2.7424781) q[2];
sx q[2];
rz(-1.9236126) q[2];
sx q[2];
rz(2.0389563) q[2];
rz(0.73978872) q[3];
sx q[3];
rz(-0.7480015) q[3];
sx q[3];
rz(1.9960777) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
