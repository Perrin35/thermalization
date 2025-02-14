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
rz(-2.8180583) q[0];
sx q[0];
rz(0.4488332) q[0];
rz(-1.1558865) q[1];
sx q[1];
rz(-1.7855676) q[1];
sx q[1];
rz(-1.1745656) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313707) q[0];
sx q[0];
rz(-2.1315984) q[0];
sx q[0];
rz(1.1602552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6392133) q[2];
sx q[2];
rz(-1.8937773) q[2];
sx q[2];
rz(-3.0634865) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60845637) q[1];
sx q[1];
rz(-1.5357596) q[1];
sx q[1];
rz(-0.94989454) q[1];
x q[2];
rz(-2.6277718) q[3];
sx q[3];
rz(-1.658545) q[3];
sx q[3];
rz(-2.3140571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.55149469) q[2];
sx q[2];
rz(-1.8415201) q[2];
sx q[2];
rz(-0.91280118) q[2];
rz(1.270795) q[3];
sx q[3];
rz(-1.0400892) q[3];
sx q[3];
rz(0.84996581) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1368644) q[0];
sx q[0];
rz(-2.5623463) q[0];
sx q[0];
rz(-2.7194523) q[0];
rz(0.36960754) q[1];
sx q[1];
rz(-2.044675) q[1];
sx q[1];
rz(0.94013989) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5053609) q[0];
sx q[0];
rz(-2.8121037) q[0];
sx q[0];
rz(-1.5785274) q[0];
rz(1.3708417) q[2];
sx q[2];
rz(-0.54206377) q[2];
sx q[2];
rz(-2.0393537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.247929) q[1];
sx q[1];
rz(-1.0808696) q[1];
sx q[1];
rz(-0.72743261) q[1];
rz(1.9056016) q[3];
sx q[3];
rz(-0.70476092) q[3];
sx q[3];
rz(-2.4089264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1993572) q[2];
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
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9966499) q[0];
sx q[0];
rz(-3.0746089) q[0];
sx q[0];
rz(1.8970733) q[0];
rz(-2.7527346) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(-0.32274524) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6038216) q[0];
sx q[0];
rz(-0.67089283) q[0];
sx q[0];
rz(0.18108271) q[0];
rz(-1.0110485) q[2];
sx q[2];
rz(-1.8343535) q[2];
sx q[2];
rz(-1.4562869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6420329) q[1];
sx q[1];
rz(-2.316058) q[1];
sx q[1];
rz(0.31742974) q[1];
rz(-pi) q[2];
rz(-1.5211608) q[3];
sx q[3];
rz(-1.2200886) q[3];
sx q[3];
rz(1.3405104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2286223) q[2];
sx q[2];
rz(-1.3801489) q[2];
sx q[2];
rz(-1.3529533) q[2];
rz(-0.27255034) q[3];
sx q[3];
rz(-1.8504146) q[3];
sx q[3];
rz(-1.4636309) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6005741) q[0];
sx q[0];
rz(-0.32222846) q[0];
sx q[0];
rz(-1.898265) q[0];
rz(2.2728032) q[1];
sx q[1];
rz(-2.181535) q[1];
sx q[1];
rz(1.0702081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58067034) q[0];
sx q[0];
rz(-1.0698478) q[0];
sx q[0];
rz(0.83606677) q[0];
rz(-pi) q[1];
rz(-1.1203499) q[2];
sx q[2];
rz(-0.82847682) q[2];
sx q[2];
rz(2.7591005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94302109) q[1];
sx q[1];
rz(-2.7987438) q[1];
sx q[1];
rz(-1.2238316) q[1];
rz(0.34262212) q[3];
sx q[3];
rz(-1.7512158) q[3];
sx q[3];
rz(0.29553771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.59226817) q[2];
sx q[2];
rz(-0.78846875) q[2];
sx q[2];
rz(0.92174021) q[2];
rz(2.8978469) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(-2.8488979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.155596) q[0];
sx q[0];
rz(-1.8159741) q[0];
sx q[0];
rz(2.6407114) q[0];
rz(-1.1174551) q[1];
sx q[1];
rz(-1.3331579) q[1];
sx q[1];
rz(0.31455988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524744) q[0];
sx q[0];
rz(-2.6388612) q[0];
sx q[0];
rz(-2.377654) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5273845) q[2];
sx q[2];
rz(-1.566326) q[2];
sx q[2];
rz(-2.4624766) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52175888) q[1];
sx q[1];
rz(-1.17621) q[1];
sx q[1];
rz(2.16741) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0231832) q[3];
sx q[3];
rz(-1.1573912) q[3];
sx q[3];
rz(-1.6550696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.306281) q[2];
sx q[2];
rz(-1.1113144) q[2];
sx q[2];
rz(-1.8446946) q[2];
rz(-0.48526192) q[3];
sx q[3];
rz(-1.5819712) q[3];
sx q[3];
rz(-1.1135134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6473963) q[0];
sx q[0];
rz(-1.2911456) q[0];
sx q[0];
rz(0.099763481) q[0];
rz(-1.8708723) q[1];
sx q[1];
rz(-0.89884177) q[1];
sx q[1];
rz(1.7521923) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4609328) q[0];
sx q[0];
rz(-1.7353199) q[0];
sx q[0];
rz(-2.8712832) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9676081) q[2];
sx q[2];
rz(-1.4672888) q[2];
sx q[2];
rz(-3.0826718) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8627117) q[1];
sx q[1];
rz(-2.2219255) q[1];
sx q[1];
rz(0.79336262) q[1];
x q[2];
rz(-1.4445969) q[3];
sx q[3];
rz(-2.3499103) q[3];
sx q[3];
rz(-2.0711957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6483267) q[2];
sx q[2];
rz(-2.6706225) q[2];
sx q[2];
rz(-0.89034447) q[2];
rz(1.9503615) q[3];
sx q[3];
rz(-1.8343364) q[3];
sx q[3];
rz(-2.2395535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4297727) q[0];
sx q[0];
rz(-3.0968102) q[0];
sx q[0];
rz(3.1088767) q[0];
rz(2.6383767) q[1];
sx q[1];
rz(-2.0423753) q[1];
sx q[1];
rz(0.12282664) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5635934) q[0];
sx q[0];
rz(-1.2482572) q[0];
sx q[0];
rz(-2.3911227) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49994464) q[2];
sx q[2];
rz(-1.2703708) q[2];
sx q[2];
rz(-2.6821313) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2715151) q[1];
sx q[1];
rz(-0.48174274) q[1];
sx q[1];
rz(-1.355257) q[1];
rz(-pi) q[2];
rz(-0.34103532) q[3];
sx q[3];
rz(-0.37217316) q[3];
sx q[3];
rz(0.91590524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7922625) q[2];
sx q[2];
rz(-0.96668875) q[2];
sx q[2];
rz(1.1518504) q[2];
rz(2.2564015) q[3];
sx q[3];
rz(-1.3693634) q[3];
sx q[3];
rz(-1.7159897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18871466) q[0];
sx q[0];
rz(-2.2080244) q[0];
sx q[0];
rz(2.2123912) q[0];
rz(-0.62905351) q[1];
sx q[1];
rz(-1.355143) q[1];
sx q[1];
rz(-2.3984875) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5261054) q[0];
sx q[0];
rz(-2.0172523) q[0];
sx q[0];
rz(3.0541712) q[0];
rz(1.6665206) q[2];
sx q[2];
rz(-2.0572955) q[2];
sx q[2];
rz(-0.41991389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.656132) q[1];
sx q[1];
rz(-0.4194057) q[1];
sx q[1];
rz(-0.16310435) q[1];
rz(-0.93502126) q[3];
sx q[3];
rz(-1.9776157) q[3];
sx q[3];
rz(-1.1233028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5762081) q[2];
sx q[2];
rz(-2.3614063) q[2];
sx q[2];
rz(2.4156477) q[2];
rz(2.7706326) q[3];
sx q[3];
rz(-1.5434664) q[3];
sx q[3];
rz(-2.7788739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1127204) q[0];
sx q[0];
rz(-2.0322587) q[0];
sx q[0];
rz(-2.0966778) q[0];
rz(-pi) q[1];
rz(-0.60301247) q[2];
sx q[2];
rz(-2.6508287) q[2];
sx q[2];
rz(-3.1163851) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5919784) q[1];
sx q[1];
rz(-1.5266298) q[1];
sx q[1];
rz(0.43818922) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56470676) q[3];
sx q[3];
rz(-2.2032524) q[3];
sx q[3];
rz(-0.45470995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5138862) q[2];
sx q[2];
rz(-2.0676925) q[2];
sx q[2];
rz(0.34454301) q[2];
rz(0.44376093) q[3];
sx q[3];
rz(-2.0659476) q[3];
sx q[3];
rz(-1.3226604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76170707) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(-2.479582) q[0];
rz(-2.8061891) q[1];
sx q[1];
rz(-1.6796203) q[1];
sx q[1];
rz(-2.2047156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59360817) q[0];
sx q[0];
rz(-2.4146705) q[0];
sx q[0];
rz(-1.0127823) q[0];
rz(1.30055) q[2];
sx q[2];
rz(-1.5480546) q[2];
sx q[2];
rz(-1.8128527) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42439991) q[1];
sx q[1];
rz(-0.48312995) q[1];
sx q[1];
rz(0.65349726) q[1];
rz(-pi) q[2];
rz(0.85031894) q[3];
sx q[3];
rz(-2.8788044) q[3];
sx q[3];
rz(0.14498728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4235437) q[2];
sx q[2];
rz(-2.9073538) q[2];
sx q[2];
rz(-2.6424109) q[2];
rz(-1.4623803) q[3];
sx q[3];
rz(-1.9951818) q[3];
sx q[3];
rz(-1.6818989) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41596169) q[0];
sx q[0];
rz(-1.9209296) q[0];
sx q[0];
rz(-0.85335535) q[0];
rz(-0.60494963) q[1];
sx q[1];
rz(-1.1440944) q[1];
sx q[1];
rz(0.49857421) q[1];
rz(-0.5729948) q[2];
sx q[2];
rz(-1.6934494) q[2];
sx q[2];
rz(2.1394503) q[2];
rz(0.50157401) q[3];
sx q[3];
rz(-1.7324823) q[3];
sx q[3];
rz(-1.7236568) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
