OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7473937) q[0];
sx q[0];
rz(-2.6497901) q[0];
sx q[0];
rz(2.9536182) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(4.6586577) q[1];
sx q[1];
rz(12.933856) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7029019) q[0];
sx q[0];
rz(-2.5950948) q[0];
sx q[0];
rz(-1.370907) q[0];
rz(-1.3221402) q[2];
sx q[2];
rz(-0.50422943) q[2];
sx q[2];
rz(0.76489514) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5012813) q[1];
sx q[1];
rz(-0.2959364) q[1];
sx q[1];
rz(0.96247767) q[1];
rz(-pi) q[2];
rz(2.8005881) q[3];
sx q[3];
rz(-1.4031646) q[3];
sx q[3];
rz(-2.117702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.964103) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.8252385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47857639) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(-2.7117803) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(0.93634161) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443003) q[0];
sx q[0];
rz(-2.8951277) q[0];
sx q[0];
rz(-0.36578567) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41994862) q[2];
sx q[2];
rz(-2.311085) q[2];
sx q[2];
rz(0.74479693) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97570005) q[1];
sx q[1];
rz(-2.4411538) q[1];
sx q[1];
rz(0.16209929) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0735047) q[3];
sx q[3];
rz(-1.5912676) q[3];
sx q[3];
rz(-2.3308144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3669746) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(-2.7152087) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(3.1085076) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8957829) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(-0.93908969) q[0];
rz(0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(-2.5476707) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1547326) q[0];
sx q[0];
rz(-1.8911456) q[0];
sx q[0];
rz(0.31517584) q[0];
rz(-0.057815876) q[2];
sx q[2];
rz(-0.88314344) q[2];
sx q[2];
rz(-2.2857894) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2908823) q[1];
sx q[1];
rz(-1.9899273) q[1];
sx q[1];
rz(-2.9078729) q[1];
rz(-1.9150312) q[3];
sx q[3];
rz(-1.1262745) q[3];
sx q[3];
rz(-1.3821186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64017355) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(1.4397941) q[2];
rz(2.7539608) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(-0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.7664292) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(2.638812) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(0.75685135) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8096136) q[0];
sx q[0];
rz(-2.6410714) q[0];
sx q[0];
rz(-2.39397) q[0];
rz(0.012979522) q[2];
sx q[2];
rz(-2.0364967) q[2];
sx q[2];
rz(-0.73060689) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0087352) q[1];
sx q[1];
rz(-0.99616226) q[1];
sx q[1];
rz(-0.43032129) q[1];
rz(-2.5960856) q[3];
sx q[3];
rz(-0.81199284) q[3];
sx q[3];
rz(0.10520392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42671529) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(-1.4871917) q[2];
rz(2.5590844) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-2.0563545) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(-1.0505189) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3462853) q[0];
sx q[0];
rz(-1.2588358) q[0];
sx q[0];
rz(-0.26766582) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1804785) q[2];
sx q[2];
rz(-2.4498307) q[2];
sx q[2];
rz(1.320653) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1861021) q[1];
sx q[1];
rz(-1.0483861) q[1];
sx q[1];
rz(-2.2239457) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0195144) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(0.36908484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22333764) q[2];
sx q[2];
rz(-2.788322) q[2];
sx q[2];
rz(-0.63344947) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(-2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(1.6075915) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(-1.4621428) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3685062) q[0];
sx q[0];
rz(-0.22531548) q[0];
sx q[0];
rz(-0.36264514) q[0];
rz(-1.2573104) q[2];
sx q[2];
rz(-0.67337155) q[2];
sx q[2];
rz(1.9999258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6219382) q[1];
sx q[1];
rz(-2.7192273) q[1];
sx q[1];
rz(-0.42610355) q[1];
rz(-pi) q[2];
rz(2.8263894) q[3];
sx q[3];
rz(-1.8990371) q[3];
sx q[3];
rz(-2.4344276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9399461) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(-0.80491006) q[2];
rz(-1.9645875) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2729623) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(-1.0246798) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(-1.0120846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597848) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(-2.3083789) q[0];
x q[1];
rz(-0.32832844) q[2];
sx q[2];
rz(-2.7331181) q[2];
sx q[2];
rz(0.35818737) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95841366) q[1];
sx q[1];
rz(-1.5378386) q[1];
sx q[1];
rz(0.54380137) q[1];
x q[2];
rz(-0.55862553) q[3];
sx q[3];
rz(-1.5010251) q[3];
sx q[3];
rz(0.39486265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53081375) q[2];
sx q[2];
rz(-1.4799708) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(-0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124509) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(-2.9274143) q[0];
rz(-2.0902436) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(-0.28373757) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9203628) q[0];
sx q[0];
rz(-1.8717248) q[0];
sx q[0];
rz(-1.5350585) q[0];
rz(-2.0869414) q[2];
sx q[2];
rz(-2.7661341) q[2];
sx q[2];
rz(1.6106538) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1188018) q[1];
sx q[1];
rz(-2.6036501) q[1];
sx q[1];
rz(-0.38938896) q[1];
rz(-pi) q[2];
rz(2.489336) q[3];
sx q[3];
rz(-2.1144923) q[3];
sx q[3];
rz(0.039747681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6909137) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(1.696375) q[2];
rz(-1.5444267) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(0.63672367) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(-3.0723363) q[0];
rz(-1.6537369) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(-1.5690631) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7135895) q[0];
sx q[0];
rz(-2.028095) q[0];
sx q[0];
rz(0.86488117) q[0];
x q[1];
rz(0.23937978) q[2];
sx q[2];
rz(-2.8068672) q[2];
sx q[2];
rz(-2.7811546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3410586) q[1];
sx q[1];
rz(-1.3375999) q[1];
sx q[1];
rz(-1.6569767) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8435555) q[3];
sx q[3];
rz(-0.66686224) q[3];
sx q[3];
rz(2.42815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1853603) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(2.774003) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-0.39127025) q[0];
sx q[0];
rz(-0.64176732) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-0.26783255) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3788911) q[0];
sx q[0];
rz(-2.2305616) q[0];
sx q[0];
rz(-2.1428109) q[0];
rz(-pi) q[1];
x q[1];
rz(3.015976) q[2];
sx q[2];
rz(-1.7526502) q[2];
sx q[2];
rz(-0.4609209) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99991998) q[1];
sx q[1];
rz(-1.281764) q[1];
sx q[1];
rz(-0.33860597) q[1];
x q[2];
rz(-0.28585163) q[3];
sx q[3];
rz(-2.0945858) q[3];
sx q[3];
rz(1.3956192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(-0.26633513) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(2.5951071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289566) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(2.3616882) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(-0.65200381) q[2];
sx q[2];
rz(-2.3163788) q[2];
sx q[2];
rz(-0.083995081) q[2];
rz(-1.7846617) q[3];
sx q[3];
rz(-0.90448096) q[3];
sx q[3];
rz(2.7703551) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
