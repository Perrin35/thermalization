OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(2.6497901) q[0];
sx q[0];
rz(9.2368035) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(-1.6245276) q[1];
sx q[1];
rz(-2.7741073) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039283218) q[0];
sx q[0];
rz(-1.4674205) q[0];
sx q[0];
rz(1.0331868) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3221402) q[2];
sx q[2];
rz(-2.6373632) q[2];
sx q[2];
rz(-0.76489514) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.130598) q[1];
sx q[1];
rz(-1.3291318) q[1];
sx q[1];
rz(0.17250891) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3931307) q[3];
sx q[3];
rz(-1.9068309) q[3];
sx q[3];
rz(2.6538268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1774896) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(-1.8356813) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.8252385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47857639) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(2.6696894) q[0];
rz(-2.7117803) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(2.205251) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443003) q[0];
sx q[0];
rz(-0.24646491) q[0];
sx q[0];
rz(0.36578567) q[0];
x q[1];
rz(-1.1510552) q[2];
sx q[2];
rz(-0.83101666) q[2];
sx q[2];
rz(-1.8120399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1658926) q[1];
sx q[1];
rz(-0.70043889) q[1];
sx q[1];
rz(0.16209929) q[1];
rz(-0.02336054) q[3];
sx q[3];
rz(-1.0682032) q[3];
sx q[3];
rz(0.7487637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77461809) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(-2.7152087) q[2];
rz(-1.2373699) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24580978) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(2.202503) q[0];
rz(-2.242873) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(0.59392196) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31375162) q[0];
sx q[0];
rz(-1.8694287) q[0];
sx q[0];
rz(-1.9065501) q[0];
rz(-1.6410286) q[2];
sx q[2];
rz(-2.4519081) q[2];
sx q[2];
rz(0.94674142) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9582639) q[1];
sx q[1];
rz(-1.7839583) q[1];
sx q[1];
rz(1.1413241) q[1];
x q[2];
rz(2.6731554) q[3];
sx q[3];
rz(-1.2611946) q[3];
sx q[3];
rz(-0.34164159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5014191) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(1.4397941) q[2];
rz(0.38763186) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(0.38813996) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(-0.50278062) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(0.75685135) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216953) q[0];
sx q[0];
rz(-1.2384402) q[0];
sx q[0];
rz(2.7601348) q[0];
x q[1];
rz(3.1286131) q[2];
sx q[2];
rz(-2.0364967) q[2];
sx q[2];
rz(-2.4109858) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.13285747) q[1];
sx q[1];
rz(-0.99616226) q[1];
sx q[1];
rz(-0.43032129) q[1];
x q[2];
rz(1.0700978) q[3];
sx q[3];
rz(-0.90161937) q[3];
sx q[3];
rz(2.5239528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42671529) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(-1.654401) q[2];
rz(2.5590844) q[3];
sx q[3];
rz(-1.094386) q[3];
sx q[3];
rz(-0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.29397598) q[0];
sx q[0];
rz(-2.0563545) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(-1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(1.0505189) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61705631) q[0];
sx q[0];
rz(-0.40818383) q[0];
sx q[0];
rz(0.88390669) q[0];
rz(-pi) q[1];
rz(2.6987223) q[2];
sx q[2];
rz(-2.1211229) q[2];
sx q[2];
rz(2.0572822) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1861021) q[1];
sx q[1];
rz(-2.0932066) q[1];
sx q[1];
rz(2.2239457) q[1];
rz(-2.3349808) q[3];
sx q[3];
rz(-1.1685373) q[3];
sx q[3];
rz(2.3276687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.918255) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(2.5081432) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(2.4244394) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69960064) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(2.8884086) q[0];
rz(-1.6075915) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(1.4621428) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77308649) q[0];
sx q[0];
rz(-0.22531548) q[0];
sx q[0];
rz(2.7789475) q[0];
x q[1];
rz(1.2573104) q[2];
sx q[2];
rz(-0.67337155) q[2];
sx q[2];
rz(1.1416669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.05789214) q[1];
sx q[1];
rz(-1.18827) q[1];
sx q[1];
rz(-1.3871357) q[1];
rz(-pi) q[2];
rz(0.83232371) q[3];
sx q[3];
rz(-0.451085) q[3];
sx q[3];
rz(3.0576599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2016466) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(0.80491006) q[2];
rz(-1.9645875) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.2729623) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(2.1169128) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(1.0120846) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1597848) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(-0.83321379) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7527533) q[2];
sx q[2];
rz(-1.4423587) q[2];
sx q[2];
rz(-1.5156137) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.55793563) q[1];
sx q[1];
rz(-0.54469889) q[1];
sx q[1];
rz(-0.063636585) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4885694) q[3];
sx q[3];
rz(-2.1279018) q[3];
sx q[3];
rz(-1.132387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.53081375) q[2];
sx q[2];
rz(-1.4799708) q[2];
sx q[2];
rz(2.7116595) q[2];
rz(-1.0144462) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124509) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(2.9274143) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(-0.28373757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8003214) q[0];
sx q[0];
rz(-2.8386136) q[0];
sx q[0];
rz(-3.0269701) q[0];
rz(-pi) q[1];
rz(-0.19212888) q[2];
sx q[2];
rz(-1.2461975) q[2];
sx q[2];
rz(1.0629551) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2092065) q[1];
sx q[1];
rz(-1.7665518) q[1];
sx q[1];
rz(-0.50435658) q[1];
x q[2];
rz(-2.3582718) q[3];
sx q[3];
rz(-0.82286994) q[3];
sx q[3];
rz(2.2058723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.45067898) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(1.696375) q[2];
rz(1.5971659) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63672367) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(0.069256393) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5690631) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7135895) q[0];
sx q[0];
rz(-2.028095) q[0];
sx q[0];
rz(-0.86488117) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8157418) q[2];
sx q[2];
rz(-1.648765) q[2];
sx q[2];
rz(-2.1577912) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8005341) q[1];
sx q[1];
rz(-1.3375999) q[1];
sx q[1];
rz(-1.6569767) q[1];
rz(2.2194355) q[3];
sx q[3];
rz(-1.7382009) q[3];
sx q[3];
rz(2.5006014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1853603) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(-2.774003) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(-0.99115133) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(1.2311252) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(0.26783255) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7627015) q[0];
sx q[0];
rz(-2.2305616) q[0];
sx q[0];
rz(-2.1428109) q[0];
x q[1];
rz(-0.12561663) q[2];
sx q[2];
rz(-1.7526502) q[2];
sx q[2];
rz(-0.4609209) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2512868) q[1];
sx q[1];
rz(-0.44154134) q[1];
sx q[1];
rz(-2.4114386) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1166499) q[3];
sx q[3];
rz(-2.5513259) q[3];
sx q[3];
rz(-2.2772307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(2.8752575) q[3];
sx q[3];
rz(-2.895152) q[3];
sx q[3];
rz(-2.5951071) q[3];
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
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(2.3616882) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(2.1521679) q[2];
sx q[2];
rz(-2.1944254) q[2];
sx q[2];
rz(2.2133322) q[2];
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
