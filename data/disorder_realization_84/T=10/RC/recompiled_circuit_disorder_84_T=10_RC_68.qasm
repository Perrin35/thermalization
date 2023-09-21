OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.26602715) q[0];
sx q[0];
rz(-0.53524435) q[0];
sx q[0];
rz(0.75403655) q[0];
rz(0.79020483) q[1];
sx q[1];
rz(-1.2269998) q[1];
sx q[1];
rz(-1.1608646) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95522674) q[0];
sx q[0];
rz(-0.47954924) q[0];
sx q[0];
rz(0.064155302) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0487469) q[2];
sx q[2];
rz(-1.4266326) q[2];
sx q[2];
rz(-2.521487) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7272211) q[1];
sx q[1];
rz(-1.1611658) q[1];
sx q[1];
rz(2.1280671) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74673064) q[3];
sx q[3];
rz(-0.24198469) q[3];
sx q[3];
rz(-0.63499588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.2966195) q[2];
sx q[2];
rz(0.66649246) q[2];
rz(-0.50981057) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(1.8681017) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78012413) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(2.0289039) q[0];
rz(2.9878222) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(-1.8033093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7897658) q[0];
sx q[0];
rz(-1.7625945) q[0];
sx q[0];
rz(-0.70621323) q[0];
x q[1];
rz(-0.43054994) q[2];
sx q[2];
rz(-0.88532788) q[2];
sx q[2];
rz(-1.4677043) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2087304) q[1];
sx q[1];
rz(-1.4534229) q[1];
sx q[1];
rz(2.7670303) q[1];
rz(-pi) q[2];
rz(-0.58740939) q[3];
sx q[3];
rz(-1.4360768) q[3];
sx q[3];
rz(0.87232529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.082836941) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(-1.178297) q[2];
rz(2.1792049) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96175471) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(1.5270365) q[0];
rz(2.4987192) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(-0.33338526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9133559) q[0];
sx q[0];
rz(-2.2524999) q[0];
sx q[0];
rz(-0.26716726) q[0];
rz(1.8893759) q[2];
sx q[2];
rz(-0.40553667) q[2];
sx q[2];
rz(0.47682724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9874939) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(0.091847329) q[1];
x q[2];
rz(-1.3473347) q[3];
sx q[3];
rz(-2.1922605) q[3];
sx q[3];
rz(0.37924757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(-2.3392759) q[2];
rz(0.24117593) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(0.10087092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872831) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(-2.5774082) q[0];
rz(2.5634649) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(-0.50813466) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80721891) q[0];
sx q[0];
rz(-2.5426572) q[0];
sx q[0];
rz(-0.81143023) q[0];
x q[1];
rz(1.5490683) q[2];
sx q[2];
rz(-2.2086775) q[2];
sx q[2];
rz(-2.4456172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9307738) q[1];
sx q[1];
rz(-2.3465119) q[1];
sx q[1];
rz(-2.7706183) q[1];
rz(-2.9594801) q[3];
sx q[3];
rz(-1.84387) q[3];
sx q[3];
rz(-1.207721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4529139) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(0.63344669) q[2];
rz(-2.5417035) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(1.6413123) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3863581) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(2.9454943) q[0];
rz(-1.8798937) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(2.0702147) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8349985) q[0];
sx q[0];
rz(-0.95933611) q[0];
sx q[0];
rz(-2.6131265) q[0];
rz(-pi) q[1];
rz(-2.7296138) q[2];
sx q[2];
rz(-2.2569071) q[2];
sx q[2];
rz(-1.0678604) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(2.1746219) q[1];
x q[2];
rz(-0.82810546) q[3];
sx q[3];
rz(-2.5525186) q[3];
sx q[3];
rz(0.38618726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1065958) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(-2.1595947) q[2];
rz(-0.18520959) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380602) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(1.8416587) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(-2.9249654) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43603555) q[0];
sx q[0];
rz(-2.2910121) q[0];
sx q[0];
rz(2.7869422) q[0];
rz(-pi) q[1];
rz(-1.8703307) q[2];
sx q[2];
rz(-0.96200633) q[2];
sx q[2];
rz(-1.2671721) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12946373) q[1];
sx q[1];
rz(-0.33540091) q[1];
sx q[1];
rz(1.5458376) q[1];
rz(0.19200237) q[3];
sx q[3];
rz(-0.90930206) q[3];
sx q[3];
rz(2.0249174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4914322) q[2];
sx q[2];
rz(-1.5101134) q[2];
sx q[2];
rz(-2.1177297) q[2];
rz(2.2653545) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(-1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.4890471) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(0.52893692) q[0];
rz(-1.5286998) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(-2.0524009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80597164) q[0];
sx q[0];
rz(-1.5651337) q[0];
sx q[0];
rz(-1.2732768) q[0];
x q[1];
rz(0.94605298) q[2];
sx q[2];
rz(-2.3240528) q[2];
sx q[2];
rz(2.3600876) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72291259) q[1];
sx q[1];
rz(-1.6039679) q[1];
sx q[1];
rz(1.564371) q[1];
x q[2];
rz(1.7004622) q[3];
sx q[3];
rz(-1.6607758) q[3];
sx q[3];
rz(2.9699096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.9160697) q[2];
rz(-1.6493753) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(-2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4847223) q[0];
sx q[0];
rz(-0.06198922) q[0];
sx q[0];
rz(2.2739676) q[0];
rz(0.067226974) q[1];
sx q[1];
rz(-2.1104689) q[1];
sx q[1];
rz(2.9464088) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0941489) q[0];
sx q[0];
rz(-0.35956811) q[0];
sx q[0];
rz(-2.5387499) q[0];
x q[1];
rz(-1.5848265) q[2];
sx q[2];
rz(-1.9719035) q[2];
sx q[2];
rz(-2.6566128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6703549) q[1];
sx q[1];
rz(-2.0002623) q[1];
sx q[1];
rz(0.045713748) q[1];
x q[2];
rz(2.04625) q[3];
sx q[3];
rz(-1.2519149) q[3];
sx q[3];
rz(1.9441324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8937257) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(2.2360738) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(0.15914966) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4147707) q[0];
sx q[0];
rz(-1.0960217) q[0];
sx q[0];
rz(1.0409521) q[0];
rz(3.0629311) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(2.7862766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61475635) q[0];
sx q[0];
rz(-0.38072452) q[0];
sx q[0];
rz(1.3044796) q[0];
rz(-pi) q[1];
rz(-1.1460733) q[2];
sx q[2];
rz(-2.6051084) q[2];
sx q[2];
rz(1.9076965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2229059) q[1];
sx q[1];
rz(-0.32864535) q[1];
sx q[1];
rz(-1.7463643) q[1];
x q[2];
rz(0.26852946) q[3];
sx q[3];
rz(-0.61754698) q[3];
sx q[3];
rz(2.7919046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(-1.7379649) q[2];
rz(-3.1353531) q[3];
sx q[3];
rz(-1.3429567) q[3];
sx q[3];
rz(-1.6041554) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733317) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(1.8433174) q[0];
rz(-2.6955993) q[1];
sx q[1];
rz(-1.1152209) q[1];
sx q[1];
rz(-0.30074063) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3215799) q[0];
sx q[0];
rz(-1.0386779) q[0];
sx q[0];
rz(-2.8768455) q[0];
rz(-pi) q[1];
rz(2.5228595) q[2];
sx q[2];
rz(-0.64090568) q[2];
sx q[2];
rz(-0.59782019) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3831331) q[1];
sx q[1];
rz(-0.29168318) q[1];
sx q[1];
rz(0.056706927) q[1];
x q[2];
rz(2.7281076) q[3];
sx q[3];
rz(-0.95260145) q[3];
sx q[3];
rz(1.2234883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62844244) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(1.1395617) q[2];
rz(1.3509753) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448485) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(-1.7383472) q[1];
sx q[1];
rz(-1.9156024) q[1];
sx q[1];
rz(1.8617873) q[1];
rz(3.1115361) q[2];
sx q[2];
rz(-3.0855784) q[2];
sx q[2];
rz(-2.2655178) q[2];
rz(-2.1915477) q[3];
sx q[3];
rz(-1.3552356) q[3];
sx q[3];
rz(-3.1380359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];