OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5300712) q[0];
sx q[0];
rz(-2.0687215) q[0];
sx q[0];
rz(-1.7142417) q[0];
rz(0.99675769) q[1];
sx q[1];
rz(-2.5280894) q[1];
sx q[1];
rz(0.065453425) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72026996) q[0];
sx q[0];
rz(-1.3822868) q[0];
sx q[0];
rz(2.6557198) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10721389) q[2];
sx q[2];
rz(-2.7709024) q[2];
sx q[2];
rz(0.057518596) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1004384) q[1];
sx q[1];
rz(-1.2693431) q[1];
sx q[1];
rz(2.235982) q[1];
x q[2];
rz(2.5259326) q[3];
sx q[3];
rz(-1.1634852) q[3];
sx q[3];
rz(3.1074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2685711) q[2];
sx q[2];
rz(-2.7950588) q[2];
sx q[2];
rz(1.2968501) q[2];
rz(2.0348564) q[3];
sx q[3];
rz(-0.96702558) q[3];
sx q[3];
rz(-1.4510252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1065555) q[0];
sx q[0];
rz(-2.4132001) q[0];
sx q[0];
rz(-1.6075217) q[0];
rz(2.2976177) q[1];
sx q[1];
rz(-1.6232099) q[1];
sx q[1];
rz(-0.27110505) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26626884) q[0];
sx q[0];
rz(-1.8401339) q[0];
sx q[0];
rz(1.6104417) q[0];
rz(1.9366802) q[2];
sx q[2];
rz(-0.26981631) q[2];
sx q[2];
rz(1.7797949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.067716777) q[1];
sx q[1];
rz(-0.42781204) q[1];
sx q[1];
rz(-2.907987) q[1];
x q[2];
rz(0.14996104) q[3];
sx q[3];
rz(-1.4984331) q[3];
sx q[3];
rz(-0.82673453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8035182) q[2];
sx q[2];
rz(-0.98793554) q[2];
sx q[2];
rz(0.79338497) q[2];
rz(-3.0986943) q[3];
sx q[3];
rz(-1.9148613) q[3];
sx q[3];
rz(-0.66810098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0637829) q[0];
sx q[0];
rz(-0.48352799) q[0];
sx q[0];
rz(-0.10312816) q[0];
rz(-2.8887796) q[1];
sx q[1];
rz(-1.7143071) q[1];
sx q[1];
rz(0.3140744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673156) q[0];
sx q[0];
rz(-2.0647103) q[0];
sx q[0];
rz(1.1252354) q[0];
rz(2.0531473) q[2];
sx q[2];
rz(-1.6214979) q[2];
sx q[2];
rz(-3.1239242) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1669877) q[1];
sx q[1];
rz(-1.0613975) q[1];
sx q[1];
rz(-2.2231119) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7400916) q[3];
sx q[3];
rz(-1.8253271) q[3];
sx q[3];
rz(-2.0930549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.126943) q[2];
sx q[2];
rz(-0.64121556) q[2];
sx q[2];
rz(-0.81652299) q[2];
rz(2.4294295) q[3];
sx q[3];
rz(-1.687259) q[3];
sx q[3];
rz(0.85675353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6461058) q[0];
sx q[0];
rz(-3.0480338) q[0];
sx q[0];
rz(-0.58864546) q[0];
rz(-2.7694287) q[1];
sx q[1];
rz(-1.2969505) q[1];
sx q[1];
rz(-2.1968496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0030768) q[0];
sx q[0];
rz(-0.32618615) q[0];
sx q[0];
rz(-1.7762) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0341357) q[2];
sx q[2];
rz(-0.54035181) q[2];
sx q[2];
rz(0.5896484) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85288799) q[1];
sx q[1];
rz(-2.1321223) q[1];
sx q[1];
rz(2.0821378) q[1];
x q[2];
rz(-0.1907986) q[3];
sx q[3];
rz(-1.9301629) q[3];
sx q[3];
rz(-0.99890733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9113691) q[2];
sx q[2];
rz(-1.5350124) q[2];
sx q[2];
rz(1.4851419) q[2];
rz(2.7459775) q[3];
sx q[3];
rz(-1.4916568) q[3];
sx q[3];
rz(-0.0609456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7570801) q[0];
sx q[0];
rz(-0.39329305) q[0];
sx q[0];
rz(1.8994037) q[0];
rz(3.0643265) q[1];
sx q[1];
rz(-1.9237513) q[1];
sx q[1];
rz(-0.088931106) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9327527) q[0];
sx q[0];
rz(-1.7249301) q[0];
sx q[0];
rz(-2.9754115) q[0];
rz(-pi) q[1];
rz(-2.5957683) q[2];
sx q[2];
rz(-1.6564545) q[2];
sx q[2];
rz(-0.31684722) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24708728) q[1];
sx q[1];
rz(-1.0997218) q[1];
sx q[1];
rz(0.31083409) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4375646) q[3];
sx q[3];
rz(-1.6682838) q[3];
sx q[3];
rz(2.9436802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5605805) q[2];
sx q[2];
rz(-1.3835013) q[2];
sx q[2];
rz(1.2233454) q[2];
rz(1.9129725) q[3];
sx q[3];
rz(-2.2556428) q[3];
sx q[3];
rz(-0.74738735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0864047) q[0];
sx q[0];
rz(-0.815027) q[0];
sx q[0];
rz(-2.1888457) q[0];
rz(1.791026) q[1];
sx q[1];
rz(-1.9074651) q[1];
sx q[1];
rz(-1.9780673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142147) q[0];
sx q[0];
rz(-2.3009536) q[0];
sx q[0];
rz(-1.7213836) q[0];
rz(-pi) q[1];
rz(0.59569401) q[2];
sx q[2];
rz(-1.3499463) q[2];
sx q[2];
rz(0.79923457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42683329) q[1];
sx q[1];
rz(-1.8140672) q[1];
sx q[1];
rz(-3.0882278) q[1];
rz(-pi) q[2];
rz(0.58168488) q[3];
sx q[3];
rz(-2.4497428) q[3];
sx q[3];
rz(-2.8210616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.20840883) q[2];
sx q[2];
rz(-1.9219425) q[2];
sx q[2];
rz(1.5737994) q[2];
rz(-0.21582223) q[3];
sx q[3];
rz(-2.5211771) q[3];
sx q[3];
rz(0.1568493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02483524) q[0];
sx q[0];
rz(-0.76501608) q[0];
sx q[0];
rz(1.5997973) q[0];
rz(0.4190017) q[1];
sx q[1];
rz(-1.1083138) q[1];
sx q[1];
rz(-1.3655183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9643819) q[0];
sx q[0];
rz(-0.96068495) q[0];
sx q[0];
rz(0.23680738) q[0];
x q[1];
rz(2.3602679) q[2];
sx q[2];
rz(-1.1622687) q[2];
sx q[2];
rz(2.9303826) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.61907459) q[1];
sx q[1];
rz(-0.72942299) q[1];
sx q[1];
rz(1.8202684) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98181866) q[3];
sx q[3];
rz(-1.2756299) q[3];
sx q[3];
rz(-1.3491614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2447723) q[2];
sx q[2];
rz(-0.12910566) q[2];
sx q[2];
rz(-1.4531892) q[2];
rz(-1.3990654) q[3];
sx q[3];
rz(-1.941967) q[3];
sx q[3];
rz(-0.80756584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.52706611) q[0];
sx q[0];
rz(-0.87166059) q[0];
sx q[0];
rz(2.6283619) q[0];
rz(-2.1649583) q[1];
sx q[1];
rz(-0.66895023) q[1];
sx q[1];
rz(-1.4248779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1126522) q[0];
sx q[0];
rz(-1.4412349) q[0];
sx q[0];
rz(0.86973637) q[0];
x q[1];
rz(3.0216359) q[2];
sx q[2];
rz(-0.55177125) q[2];
sx q[2];
rz(-0.063869501) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8911881) q[1];
sx q[1];
rz(-2.0051167) q[1];
sx q[1];
rz(-2.8211947) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1838328) q[3];
sx q[3];
rz(-2.9385205) q[3];
sx q[3];
rz(1.9436294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49154115) q[2];
sx q[2];
rz(-1.0074793) q[2];
sx q[2];
rz(-0.36637351) q[2];
rz(0.63792396) q[3];
sx q[3];
rz(-1.7584691) q[3];
sx q[3];
rz(-1.6320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76846182) q[0];
sx q[0];
rz(-2.8167384) q[0];
sx q[0];
rz(-1.1573867) q[0];
rz(-0.53818446) q[1];
sx q[1];
rz(-1.8073852) q[1];
sx q[1];
rz(-0.023712637) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2185604) q[0];
sx q[0];
rz(-2.5171748) q[0];
sx q[0];
rz(3.0397281) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3228119) q[2];
sx q[2];
rz(-1.2597558) q[2];
sx q[2];
rz(-3.0672665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38172028) q[1];
sx q[1];
rz(-0.64331912) q[1];
sx q[1];
rz(2.5030959) q[1];
x q[2];
rz(1.0162418) q[3];
sx q[3];
rz(-0.23532495) q[3];
sx q[3];
rz(-0.91431844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9862765) q[2];
sx q[2];
rz(-1.3331058) q[2];
sx q[2];
rz(-0.39815608) q[2];
rz(2.6768173) q[3];
sx q[3];
rz(-2.0538752) q[3];
sx q[3];
rz(-1.3409415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.583113) q[0];
sx q[0];
rz(-2.8093503) q[0];
sx q[0];
rz(0.94859052) q[0];
rz(0.96725431) q[1];
sx q[1];
rz(-1.1525258) q[1];
sx q[1];
rz(2.1327877) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0111158) q[0];
sx q[0];
rz(-1.563894) q[0];
sx q[0];
rz(2.3409977) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1494184) q[2];
sx q[2];
rz(-2.1632901) q[2];
sx q[2];
rz(-0.078607056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12647835) q[1];
sx q[1];
rz(-0.51869828) q[1];
sx q[1];
rz(-0.37094231) q[1];
rz(-1.7901374) q[3];
sx q[3];
rz(-2.2171341) q[3];
sx q[3];
rz(-1.5462033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7271416) q[2];
sx q[2];
rz(-0.87697566) q[2];
sx q[2];
rz(-3.1205175) q[2];
rz(-0.4840788) q[3];
sx q[3];
rz(-1.1295854) q[3];
sx q[3];
rz(2.4626203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2912343) q[0];
sx q[0];
rz(-1.5342916) q[0];
sx q[0];
rz(1.7271484) q[0];
rz(0.5961295) q[1];
sx q[1];
rz(-2.3565751) q[1];
sx q[1];
rz(-0.48859488) q[1];
rz(-0.1487371) q[2];
sx q[2];
rz(-1.3631459) q[2];
sx q[2];
rz(2.7447328) q[2];
rz(0.10808839) q[3];
sx q[3];
rz(-1.3664403) q[3];
sx q[3];
rz(-1.4667778) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
