OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(-0.35020819) q[0];
sx q[0];
rz(-0.36663088) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(1.4555414) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35858425) q[0];
sx q[0];
rz(-1.9134247) q[0];
sx q[0];
rz(1.1230099) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86906616) q[2];
sx q[2];
rz(-1.367374) q[2];
sx q[2];
rz(0.8806526) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75016025) q[1];
sx q[1];
rz(-0.62392985) q[1];
sx q[1];
rz(-1.0990012) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5600169) q[3];
sx q[3];
rz(-1.5228776) q[3];
sx q[3];
rz(0.77424327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(2.74995) q[2];
rz(-0.13970217) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166819) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(2.2712767) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.93719) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8613974) q[0];
sx q[0];
rz(-1.9635927) q[0];
sx q[0];
rz(-1.1108206) q[0];
rz(-pi) q[1];
rz(2.9559829) q[2];
sx q[2];
rz(-0.87366548) q[2];
sx q[2];
rz(-2.6999161) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.054489) q[1];
sx q[1];
rz(-0.37885715) q[1];
sx q[1];
rz(0.27146753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8562206) q[3];
sx q[3];
rz(-2.223569) q[3];
sx q[3];
rz(-0.6245581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5045972) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(1.2191999) q[2];
rz(-2.7870264) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5511659) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(-0.61808008) q[0];
rz(0.016050054) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(1.191167) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1026099) q[0];
sx q[0];
rz(-1.257574) q[0];
sx q[0];
rz(0.15977504) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3267924) q[2];
sx q[2];
rz(-1.2952431) q[2];
sx q[2];
rz(-0.22305605) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.788113) q[1];
sx q[1];
rz(-2.8444926) q[1];
sx q[1];
rz(-0.77570813) q[1];
rz(-pi) q[2];
rz(2.5997945) q[3];
sx q[3];
rz(-1.2216611) q[3];
sx q[3];
rz(-0.19402129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1742192) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(-1.013914) q[2];
rz(-1.3570471) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8106666) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(0.20430918) q[0];
rz(-1.3775685) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(2.2185982) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92813338) q[0];
sx q[0];
rz(-1.8753578) q[0];
sx q[0];
rz(-0.34780985) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26206215) q[2];
sx q[2];
rz(-0.045496551) q[2];
sx q[2];
rz(-2.2646575) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1624565) q[1];
sx q[1];
rz(-0.41185954) q[1];
sx q[1];
rz(-1.3084175) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8664341) q[3];
sx q[3];
rz(-1.5878864) q[3];
sx q[3];
rz(-2.7279502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(-2.6320809) q[2];
rz(2.4984958) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(-2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5287857) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(0.91941419) q[0];
rz(0.98584229) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(1.7808419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0589941) q[0];
sx q[0];
rz(-1.3816125) q[0];
sx q[0];
rz(-2.9181705) q[0];
rz(-pi) q[1];
rz(-2.1941357) q[2];
sx q[2];
rz(-2.2286378) q[2];
sx q[2];
rz(0.62589494) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9779224) q[1];
sx q[1];
rz(-1.1174035) q[1];
sx q[1];
rz(1.3400673) q[1];
rz(-pi) q[2];
rz(-0.97134437) q[3];
sx q[3];
rz(-1.6347307) q[3];
sx q[3];
rz(1.5980699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(-0.11486593) q[2];
rz(-2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546656) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(-1.2448357) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(-2.8663666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6557065) q[0];
sx q[0];
rz(-2.5285892) q[0];
sx q[0];
rz(2.1711151) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.198344) q[2];
sx q[2];
rz(-2.3970282) q[2];
sx q[2];
rz(-0.17472357) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5425307) q[1];
sx q[1];
rz(-1.9179357) q[1];
sx q[1];
rz(-0.51862006) q[1];
rz(-pi) q[2];
rz(-1.4554548) q[3];
sx q[3];
rz(-2.5177258) q[3];
sx q[3];
rz(1.8231525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.747921) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(2.4300872) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(-0.59035629) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8608619) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(1.77805) q[0];
rz(1.4200312) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-2.4218959) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5170383) q[0];
sx q[0];
rz(-1.8554243) q[0];
sx q[0];
rz(2.4863003) q[0];
x q[1];
rz(2.4160956) q[2];
sx q[2];
rz(-0.91425397) q[2];
sx q[2];
rz(-0.072349116) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4622725) q[1];
sx q[1];
rz(-1.6151036) q[1];
sx q[1];
rz(-0.54507749) q[1];
rz(-2.4607055) q[3];
sx q[3];
rz(-1.3586992) q[3];
sx q[3];
rz(-1.6915481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65934962) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-0.68816319) q[2];
rz(-0.041953772) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(-0.18653175) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(-1.4950745) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625531) q[0];
sx q[0];
rz(-0.4058668) q[0];
sx q[0];
rz(2.5748475) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.375884) q[2];
sx q[2];
rz(-2.2576828) q[2];
sx q[2];
rz(1.9538823) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9547792) q[1];
sx q[1];
rz(-2.8615132) q[1];
sx q[1];
rz(-0.38444744) q[1];
x q[2];
rz(-2.5008194) q[3];
sx q[3];
rz(-0.91859222) q[3];
sx q[3];
rz(-2.3280502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24215332) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(3.126826) q[2];
rz(-1.8642558) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(-1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(0.87798464) q[0];
rz(0.19628482) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(1.7810129) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1454865) q[0];
sx q[0];
rz(-1.6084451) q[0];
sx q[0];
rz(-0.99897511) q[0];
rz(-0.94366818) q[2];
sx q[2];
rz(-1.6086279) q[2];
sx q[2];
rz(2.949122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66879241) q[1];
sx q[1];
rz(-1.2392514) q[1];
sx q[1];
rz(1.214295) q[1];
x q[2];
rz(-1.1514879) q[3];
sx q[3];
rz(-1.4290223) q[3];
sx q[3];
rz(-1.3594974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62375162) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(-0.8927792) q[2];
rz(3.1023074) q[3];
sx q[3];
rz(-1.4792484) q[3];
sx q[3];
rz(-0.75004309) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2980625) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(3.0950586) q[0];
rz(-2.2325366) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-2.4216901) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68144804) q[0];
sx q[0];
rz(-1.8143192) q[0];
sx q[0];
rz(-0.081865099) q[0];
rz(-pi) q[1];
rz(1.5762395) q[2];
sx q[2];
rz(-1.8796225) q[2];
sx q[2];
rz(-2.960161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3958017) q[1];
sx q[1];
rz(-1.5553286) q[1];
sx q[1];
rz(-1.8176816) q[1];
rz(-pi) q[2];
rz(0.7474483) q[3];
sx q[3];
rz(-2.2564853) q[3];
sx q[3];
rz(1.7789343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7211192) q[2];
sx q[2];
rz(-1.7420008) q[2];
sx q[2];
rz(0.026467888) q[2];
rz(-1.3039533) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(-1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883041) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(0.2164671) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(0.36454501) q[2];
sx q[2];
rz(-0.74082965) q[2];
sx q[2];
rz(-0.34168591) q[2];
rz(1.5218232) q[3];
sx q[3];
rz(-2.5544142) q[3];
sx q[3];
rz(2.1094473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
