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
rz(1.264313) q[0];
sx q[0];
rz(6.6067196) q[0];
sx q[0];
rz(8.9759448) q[0];
rz(-1.1558865) q[1];
sx q[1];
rz(-1.7855676) q[1];
sx q[1];
rz(-1.1745656) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1313707) q[0];
sx q[0];
rz(-1.0099942) q[0];
sx q[0];
rz(1.1602552) q[0];
rz(-pi) q[1];
rz(-0.20148142) q[2];
sx q[2];
rz(-0.32989943) q[2];
sx q[2];
rz(3.0070674) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0112746) q[1];
sx q[1];
rz(-2.5198333) q[1];
sx q[1];
rz(1.6309726) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51382084) q[3];
sx q[3];
rz(-1.4830477) q[3];
sx q[3];
rz(-0.82753554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.55149469) q[2];
sx q[2];
rz(-1.8415201) q[2];
sx q[2];
rz(0.91280118) q[2];
rz(1.8707976) q[3];
sx q[3];
rz(-1.0400892) q[3];
sx q[3];
rz(-0.84996581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1368644) q[0];
sx q[0];
rz(-2.5623463) q[0];
sx q[0];
rz(0.42214033) q[0];
rz(2.7719851) q[1];
sx q[1];
rz(-2.044675) q[1];
sx q[1];
rz(2.2014528) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62806118) q[0];
sx q[0];
rz(-1.9002751) q[0];
sx q[0];
rz(-0.0026436289) q[0];
rz(-pi) q[1];
rz(1.0375848) q[2];
sx q[2];
rz(-1.6734481) q[2];
sx q[2];
rz(-2.8449322) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9782422) q[1];
sx q[1];
rz(-0.85127318) q[1];
sx q[1];
rz(-0.6759599) q[1];
x q[2];
rz(-2.8690954) q[3];
sx q[3];
rz(-2.2292308) q[3];
sx q[3];
rz(0.30425018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9422354) q[2];
sx q[2];
rz(-2.4713559) q[2];
sx q[2];
rz(2.2352236) q[2];
rz(2.1622315) q[3];
sx q[3];
rz(-2.7370079) q[3];
sx q[3];
rz(-1.8766859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9966499) q[0];
sx q[0];
rz(-3.0746089) q[0];
sx q[0];
rz(-1.8970733) q[0];
rz(-0.38885802) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(-2.8188474) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9661315) q[0];
sx q[0];
rz(-1.6829938) q[0];
sx q[0];
rz(2.4787122) q[0];
x q[1];
rz(2.8333146) q[2];
sx q[2];
rz(-2.1090504) q[2];
sx q[2];
rz(0.27632144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14798627) q[1];
sx q[1];
rz(-1.3393511) q[1];
sx q[1];
rz(-0.79995059) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6204319) q[3];
sx q[3];
rz(-1.921504) q[3];
sx q[3];
rz(-1.3405104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2286223) q[2];
sx q[2];
rz(-1.7614438) q[2];
sx q[2];
rz(1.3529533) q[2];
rz(-2.8690423) q[3];
sx q[3];
rz(-1.291178) q[3];
sx q[3];
rz(-1.4636309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54101855) q[0];
sx q[0];
rz(-2.8193642) q[0];
sx q[0];
rz(1.2433276) q[0];
rz(2.2728032) q[1];
sx q[1];
rz(-0.96005762) q[1];
sx q[1];
rz(-1.0702081) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3995099) q[0];
sx q[0];
rz(-2.1993981) q[0];
sx q[0];
rz(-2.5058772) q[0];
rz(-pi) q[1];
rz(2.3467873) q[2];
sx q[2];
rz(-1.2441976) q[2];
sx q[2];
rz(-2.2692533) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29954545) q[1];
sx q[1];
rz(-1.6853602) q[1];
sx q[1];
rz(-1.8946527) q[1];
x q[2];
rz(0.49740187) q[3];
sx q[3];
rz(-0.38555749) q[3];
sx q[3];
rz(1.7413643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59226817) q[2];
sx q[2];
rz(-2.3531239) q[2];
sx q[2];
rz(2.2198524) q[2];
rz(-2.8978469) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(2.8488979) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.155596) q[0];
sx q[0];
rz(-1.8159741) q[0];
sx q[0];
rz(2.6407114) q[0];
rz(2.0241375) q[1];
sx q[1];
rz(-1.8084348) q[1];
sx q[1];
rz(-0.31455988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68334333) q[0];
sx q[0];
rz(-1.2309845) q[0];
sx q[0];
rz(-2.763624) q[0];
x q[1];
rz(-1.4681513) q[2];
sx q[2];
rz(-0.043641239) q[2];
sx q[2];
rz(0.99422821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.79364538) q[1];
sx q[1];
rz(-1.0255019) q[1];
sx q[1];
rz(-0.46635638) q[1];
x q[2];
rz(-3.0231832) q[3];
sx q[3];
rz(-1.9842015) q[3];
sx q[3];
rz(-1.6550696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.306281) q[2];
sx q[2];
rz(-2.0302782) q[2];
sx q[2];
rz(-1.8446946) q[2];
rz(-0.48526192) q[3];
sx q[3];
rz(-1.5596215) q[3];
sx q[3];
rz(1.1135134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6473963) q[0];
sx q[0];
rz(-1.2911456) q[0];
sx q[0];
rz(3.0418292) q[0];
rz(-1.8708723) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(-1.7521923) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4978964) q[0];
sx q[0];
rz(-0.31539105) q[0];
sx q[0];
rz(-0.55625486) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17398457) q[2];
sx q[2];
rz(-1.6743039) q[2];
sx q[2];
rz(-3.0826718) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.278881) q[1];
sx q[1];
rz(-2.2219255) q[1];
sx q[1];
rz(-0.79336262) q[1];
rz(1.4445969) q[3];
sx q[3];
rz(-2.3499103) q[3];
sx q[3];
rz(-1.070397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49326593) q[2];
sx q[2];
rz(-2.6706225) q[2];
sx q[2];
rz(2.2512482) q[2];
rz(1.1912311) q[3];
sx q[3];
rz(-1.3072562) q[3];
sx q[3];
rz(0.90203917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4297727) q[0];
sx q[0];
rz(-3.0968102) q[0];
sx q[0];
rz(-0.032715948) q[0];
rz(-0.50321594) q[1];
sx q[1];
rz(-2.0423753) q[1];
sx q[1];
rz(-3.018766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8211201) q[0];
sx q[0];
rz(-2.3373465) q[0];
sx q[0];
rz(2.6859317) q[0];
rz(-1.2314447) q[2];
sx q[2];
rz(-1.095158) q[2];
sx q[2];
rz(-0.95107691) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2715151) q[1];
sx q[1];
rz(-0.48174274) q[1];
sx q[1];
rz(1.355257) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7006247) q[3];
sx q[3];
rz(-1.2210088) q[3];
sx q[3];
rz(-1.2799124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34933019) q[2];
sx q[2];
rz(-0.96668875) q[2];
sx q[2];
rz(-1.9897423) q[2];
rz(0.88519111) q[3];
sx q[3];
rz(-1.7722292) q[3];
sx q[3];
rz(-1.7159897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952878) q[0];
sx q[0];
rz(-2.2080244) q[0];
sx q[0];
rz(-2.2123912) q[0];
rz(-0.62905351) q[1];
sx q[1];
rz(-1.355143) q[1];
sx q[1];
rz(-2.3984875) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5261054) q[0];
sx q[0];
rz(-2.0172523) q[0];
sx q[0];
rz(-3.0541712) q[0];
rz(-pi) q[1];
rz(-2.6531947) q[2];
sx q[2];
rz(-1.4862068) q[2];
sx q[2];
rz(-1.9458488) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.23452246) q[1];
sx q[1];
rz(-1.6369695) q[1];
sx q[1];
rz(-2.727134) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.943103) q[3];
sx q[3];
rz(-2.4022958) q[3];
sx q[3];
rz(-0.044667808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5653845) q[2];
sx q[2];
rz(-0.78018633) q[2];
sx q[2];
rz(-0.72594491) q[2];
rz(-0.37096008) q[3];
sx q[3];
rz(-1.5981263) q[3];
sx q[3];
rz(2.7788739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2711656) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(0.56394947) q[0];
rz(1.14934) q[1];
sx q[1];
rz(-0.96200689) q[1];
sx q[1];
rz(-2.3990778) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0288723) q[0];
sx q[0];
rz(-1.109334) q[0];
sx q[0];
rz(-1.0449148) q[0];
rz(-pi) q[1];
rz(-1.276539) q[2];
sx q[2];
rz(-1.1721436) q[2];
sx q[2];
rz(-2.5039303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5919784) q[1];
sx q[1];
rz(-1.6149628) q[1];
sx q[1];
rz(0.43818922) q[1];
rz(-pi) q[2];
rz(-2.5768859) q[3];
sx q[3];
rz(-2.2032524) q[3];
sx q[3];
rz(-2.6868827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5138862) q[2];
sx q[2];
rz(-2.0676925) q[2];
sx q[2];
rz(2.7970496) q[2];
rz(-0.44376093) q[3];
sx q[3];
rz(-2.0659476) q[3];
sx q[3];
rz(-1.8189323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3798856) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(-0.66201061) q[0];
rz(2.8061891) q[1];
sx q[1];
rz(-1.4619724) q[1];
sx q[1];
rz(0.9368771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59360817) q[0];
sx q[0];
rz(-0.72692211) q[0];
sx q[0];
rz(-1.0127823) q[0];
x q[1];
rz(1.30055) q[2];
sx q[2];
rz(-1.5480546) q[2];
sx q[2];
rz(1.3287399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8530218) q[1];
sx q[1];
rz(-1.9485546) q[1];
sx q[1];
rz(-1.8795345) q[1];
x q[2];
rz(1.7702661) q[3];
sx q[3];
rz(-1.7430309) q[3];
sx q[3];
rz(2.1290092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4235437) q[2];
sx q[2];
rz(-2.9073538) q[2];
sx q[2];
rz(0.49918175) q[2];
rz(-1.6792123) q[3];
sx q[3];
rz(-1.1464109) q[3];
sx q[3];
rz(1.4596938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.725631) q[0];
sx q[0];
rz(-1.9209296) q[0];
sx q[0];
rz(-0.85335535) q[0];
rz(-0.60494963) q[1];
sx q[1];
rz(-1.1440944) q[1];
sx q[1];
rz(0.49857421) q[1];
rz(-2.9180183) q[2];
sx q[2];
rz(-0.58453544) q[2];
sx q[2];
rz(0.75605308) q[2];
rz(-1.3868757) q[3];
sx q[3];
rz(-1.0763604) q[3];
sx q[3];
rz(3.0767783) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
