OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(6.0072748) q[0];
sx q[0];
rz(10.732565) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(-1.5712665) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64496541) q[0];
sx q[0];
rz(-2.6500406) q[0];
sx q[0];
rz(-0.77792032) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70648944) q[2];
sx q[2];
rz(-2.2405365) q[2];
sx q[2];
rz(-2.0073839) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.03304122) q[1];
sx q[1];
rz(-1.3290977) q[1];
sx q[1];
rz(2.7944195) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0516112) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(-0.68457505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87542614) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(-2.0092633) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.9448626) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(2.9247608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8502055) q[0];
sx q[0];
rz(-2.4225525) q[0];
sx q[0];
rz(-2.015381) q[0];
rz(-pi) q[1];
rz(0.88044135) q[2];
sx q[2];
rz(-2.464622) q[2];
sx q[2];
rz(-2.0073236) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.617802) q[1];
sx q[1];
rz(-0.70980598) q[1];
sx q[1];
rz(-1.2426504) q[1];
rz(0.87644491) q[3];
sx q[3];
rz(-1.0507686) q[3];
sx q[3];
rz(2.2976573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.310114) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(-2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6771616) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(0.60423869) q[0];
rz(-1.8151981) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(0.93260971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7639097) q[0];
sx q[0];
rz(-1.253486) q[0];
sx q[0];
rz(-0.056563932) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3435752) q[2];
sx q[2];
rz(-1.3933239) q[2];
sx q[2];
rz(0.15572671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.052913594) q[1];
sx q[1];
rz(-1.0148078) q[1];
sx q[1];
rz(-1.710379) q[1];
rz(2.695735) q[3];
sx q[3];
rz(-1.0599531) q[3];
sx q[3];
rz(1.0872935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3595235) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(0.50022593) q[0];
rz(0.80530986) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(-1.4979699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26101199) q[0];
sx q[0];
rz(-1.2769165) q[0];
sx q[0];
rz(-2.0902993) q[0];
x q[1];
rz(1.8565606) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(-2.6464268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8151617) q[1];
sx q[1];
rz(-1.3109428) q[1];
sx q[1];
rz(-1.8111147) q[1];
rz(2.0235396) q[3];
sx q[3];
rz(-1.9948043) q[3];
sx q[3];
rz(0.31183576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74636373) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(-2.4397819) q[2];
rz(2.3102405) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(-0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.9005301) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(0.81533122) q[0];
rz(1.6197846) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-2.0933847) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5883023) q[0];
sx q[0];
rz(-0.36670812) q[0];
sx q[0];
rz(0.81272965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65152119) q[2];
sx q[2];
rz(-1.5805575) q[2];
sx q[2];
rz(-0.62192813) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8909) q[1];
sx q[1];
rz(-2.7437468) q[1];
sx q[1];
rz(-2.489151) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0190373) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(-2.2321731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(2.0416416) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0681756) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(-0.90240479) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(-3.0117603) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5141402) q[0];
sx q[0];
rz(-2.1462626) q[0];
sx q[0];
rz(2.0155725) q[0];
rz(-pi) q[1];
rz(-2.7733299) q[2];
sx q[2];
rz(-2.1149181) q[2];
sx q[2];
rz(0.95168176) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9784769) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(2.5549868) q[1];
rz(-pi) q[2];
x q[2];
rz(1.767166) q[3];
sx q[3];
rz(-2.11103) q[3];
sx q[3];
rz(1.9062717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3123902) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(-2.9373346) q[2];
rz(-1.2060818) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7234574) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.6947421) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(-2.4553305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8780554) q[0];
sx q[0];
rz(-2.005308) q[0];
sx q[0];
rz(2.6146019) q[0];
x q[1];
rz(2.1967728) q[2];
sx q[2];
rz(-2.2959024) q[2];
sx q[2];
rz(-0.041989728) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2346748) q[1];
sx q[1];
rz(-1.7559116) q[1];
sx q[1];
rz(-0.075637416) q[1];
rz(-pi) q[2];
rz(-2.5043082) q[3];
sx q[3];
rz(-2.6316959) q[3];
sx q[3];
rz(-2.2989458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(0.0017722842) q[2];
rz(-2.5799675) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(-1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5381662) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.6400281) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9672464) q[0];
sx q[0];
rz(-1.5513199) q[0];
sx q[0];
rz(-2.0635701) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.026560606) q[2];
sx q[2];
rz(-1.5090669) q[2];
sx q[2];
rz(-0.82690566) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.20128076) q[1];
sx q[1];
rz(-0.33946013) q[1];
sx q[1];
rz(2.8708354) q[1];
rz(-pi) q[2];
rz(-0.6286962) q[3];
sx q[3];
rz(-1.0843715) q[3];
sx q[3];
rz(-1.6823671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7897196) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(1.8224576) q[2];
rz(-1.2119279) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(-0.31931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-0.35167545) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68181224) q[0];
sx q[0];
rz(-2.5373055) q[0];
sx q[0];
rz(2.2629645) q[0];
rz(-0.69182379) q[2];
sx q[2];
rz(-1.8080538) q[2];
sx q[2];
rz(-2.0123864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-2.3304686) q[1];
sx q[1];
rz(-3.0685436) q[1];
rz(1.5185235) q[3];
sx q[3];
rz(-1.4613323) q[3];
sx q[3];
rz(1.048552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.3433156) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(1.2822255) q[2];
rz(-1.4964237) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(-2.9472651) q[0];
rz(-2.1037897) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(-2.1077572) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949085) q[0];
sx q[0];
rz(-2.2241728) q[0];
sx q[0];
rz(2.9772467) q[0];
x q[1];
rz(-0.50844426) q[2];
sx q[2];
rz(-2.2061081) q[2];
sx q[2];
rz(-0.19257643) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2227877) q[1];
sx q[1];
rz(-2.5606887) q[1];
sx q[1];
rz(-1.3851628) q[1];
rz(-pi) q[2];
rz(0.25063534) q[3];
sx q[3];
rz(-2.4126629) q[3];
sx q[3];
rz(0.38754101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(2.5058084) q[2];
rz(0.27030269) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.6939659) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(1.7059965) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(2.1736017) q[2];
sx q[2];
rz(-1.597076) q[2];
sx q[2];
rz(1.1521641) q[2];
rz(-0.070449645) q[3];
sx q[3];
rz(-1.2415213) q[3];
sx q[3];
rz(0.51125676) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
