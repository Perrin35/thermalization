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
rz(2.1228696) q[0];
sx q[0];
rz(-2.2824204) q[0];
sx q[0];
rz(0.81508842) q[0];
rz(1.9563142) q[1];
sx q[1];
rz(-1.7307245) q[1];
sx q[1];
rz(-1.0676395) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0420494) q[0];
sx q[0];
rz(-1.1394412) q[0];
sx q[0];
rz(0.026148034) q[0];
rz(2.5065866) q[2];
sx q[2];
rz(-2.5189812) q[2];
sx q[2];
rz(-0.16281637) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9166482) q[1];
sx q[1];
rz(-1.4743346) q[1];
sx q[1];
rz(-1.2526082) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2514011) q[3];
sx q[3];
rz(-0.71139627) q[3];
sx q[3];
rz(-2.9044202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4472569) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(-1.2747964) q[2];
rz(2.6796807) q[3];
sx q[3];
rz(-0.67449823) q[3];
sx q[3];
rz(2.0089202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79134113) q[0];
sx q[0];
rz(-0.30838648) q[0];
sx q[0];
rz(1.8885008) q[0];
rz(0.14532267) q[1];
sx q[1];
rz(-1.7456313) q[1];
sx q[1];
rz(-1.0911509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6515132) q[0];
sx q[0];
rz(-1.7848178) q[0];
sx q[0];
rz(1.1351311) q[0];
x q[1];
rz(2.4895489) q[2];
sx q[2];
rz(-0.55527675) q[2];
sx q[2];
rz(2.9342143) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.021046358) q[1];
sx q[1];
rz(-1.7851549) q[1];
sx q[1];
rz(1.7477186) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49511045) q[3];
sx q[3];
rz(-1.2155327) q[3];
sx q[3];
rz(-2.2877467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6609409) q[2];
sx q[2];
rz(-1.3903684) q[2];
sx q[2];
rz(-0.52978984) q[2];
rz(2.3482813) q[3];
sx q[3];
rz(-1.6042234) q[3];
sx q[3];
rz(-0.62354273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6563501) q[0];
sx q[0];
rz(-0.95004496) q[0];
sx q[0];
rz(2.6128838) q[0];
rz(-0.54620019) q[1];
sx q[1];
rz(-2.1827953) q[1];
sx q[1];
rz(-0.34034696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0673235) q[0];
sx q[0];
rz(-2.7087697) q[0];
sx q[0];
rz(2.4433663) q[0];
rz(-0.23480798) q[2];
sx q[2];
rz(-1.4291414) q[2];
sx q[2];
rz(2.0675142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9593411) q[1];
sx q[1];
rz(-1.6653582) q[1];
sx q[1];
rz(-1.5675117) q[1];
x q[2];
rz(1.163655) q[3];
sx q[3];
rz(-0.7837067) q[3];
sx q[3];
rz(-0.78898417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9693552) q[2];
sx q[2];
rz(-0.41219553) q[2];
sx q[2];
rz(-2.7117512) q[2];
rz(-2.3853081) q[3];
sx q[3];
rz(-2.9679306) q[3];
sx q[3];
rz(2.3156796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38465685) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(1.6480308) q[0];
rz(-2.3736296) q[1];
sx q[1];
rz(-0.47052828) q[1];
sx q[1];
rz(-0.54642645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38619216) q[0];
sx q[0];
rz(-1.6419171) q[0];
sx q[0];
rz(-3.0767308) q[0];
rz(2.0794608) q[2];
sx q[2];
rz(-2.0360247) q[2];
sx q[2];
rz(2.0874799) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5560234) q[1];
sx q[1];
rz(-1.4735735) q[1];
sx q[1];
rz(2.9929964) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91937842) q[3];
sx q[3];
rz(-2.1438144) q[3];
sx q[3];
rz(-1.023055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7118608) q[2];
sx q[2];
rz(-1.6542566) q[2];
sx q[2];
rz(-0.30409733) q[2];
rz(-1.3652623) q[3];
sx q[3];
rz(-1.920776) q[3];
sx q[3];
rz(-2.064866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6292608) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(3.1410134) q[0];
rz(-2.0594788) q[1];
sx q[1];
rz(-0.52434701) q[1];
sx q[1];
rz(-2.2023315) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88193306) q[0];
sx q[0];
rz(-1.1355917) q[0];
sx q[0];
rz(2.4324904) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0505191) q[2];
sx q[2];
rz(-2.4726598) q[2];
sx q[2];
rz(2.875653) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56783695) q[1];
sx q[1];
rz(-1.3569965) q[1];
sx q[1];
rz(2.1478199) q[1];
rz(-pi) q[2];
rz(2.9143798) q[3];
sx q[3];
rz(-1.2649396) q[3];
sx q[3];
rz(-1.3671041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.038736343) q[2];
sx q[2];
rz(-1.2612217) q[2];
sx q[2];
rz(0.2612513) q[2];
rz(-2.2767565) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(-2.849546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84457266) q[0];
sx q[0];
rz(-1.7140056) q[0];
sx q[0];
rz(-0.20508668) q[0];
rz(2.0948441) q[1];
sx q[1];
rz(-1.3958684) q[1];
sx q[1];
rz(-0.37839016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32961938) q[0];
sx q[0];
rz(-1.9549592) q[0];
sx q[0];
rz(-0.20150082) q[0];
rz(-pi) q[1];
rz(-0.47176265) q[2];
sx q[2];
rz(-1.7182351) q[2];
sx q[2];
rz(2.9094537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4503895) q[1];
sx q[1];
rz(-0.65237633) q[1];
sx q[1];
rz(-0.97263402) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50726733) q[3];
sx q[3];
rz(-2.3823822) q[3];
sx q[3];
rz(1.4634446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4947074) q[2];
sx q[2];
rz(-1.9834221) q[2];
sx q[2];
rz(2.0873439) q[2];
rz(2.4833637) q[3];
sx q[3];
rz(-2.4963278) q[3];
sx q[3];
rz(0.48404199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14195104) q[0];
sx q[0];
rz(-1.208409) q[0];
sx q[0];
rz(-0.64055881) q[0];
rz(-2.7236252) q[1];
sx q[1];
rz(-1.2203981) q[1];
sx q[1];
rz(-2.3366065) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1520839) q[0];
sx q[0];
rz(-1.594127) q[0];
sx q[0];
rz(-0.70744608) q[0];
x q[1];
rz(-1.1518794) q[2];
sx q[2];
rz(-2.3275314) q[2];
sx q[2];
rz(-0.69413041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1272443) q[1];
sx q[1];
rz(-1.5509203) q[1];
sx q[1];
rz(-1.5991421) q[1];
rz(-pi) q[2];
rz(-2.5500357) q[3];
sx q[3];
rz(-1.8186967) q[3];
sx q[3];
rz(1.7660559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54887041) q[2];
sx q[2];
rz(-0.9950811) q[2];
sx q[2];
rz(-2.3804046) q[2];
rz(2.393764) q[3];
sx q[3];
rz(-1.8239832) q[3];
sx q[3];
rz(-0.67659155) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51650301) q[0];
sx q[0];
rz(-0.28446063) q[0];
sx q[0];
rz(-1.6647343) q[0];
rz(2.6853216) q[1];
sx q[1];
rz(-1.7218593) q[1];
sx q[1];
rz(0.87108535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3960311) q[0];
sx q[0];
rz(-2.1848688) q[0];
sx q[0];
rz(-1.5860193) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7256868) q[2];
sx q[2];
rz(-1.2300856) q[2];
sx q[2];
rz(1.0122077) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2312517) q[1];
sx q[1];
rz(-1.3518855) q[1];
sx q[1];
rz(0.33556767) q[1];
rz(1.4389478) q[3];
sx q[3];
rz(-0.42964298) q[3];
sx q[3];
rz(2.1983918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90290922) q[2];
sx q[2];
rz(-0.52672714) q[2];
sx q[2];
rz(2.5229559) q[2];
rz(3.1213308) q[3];
sx q[3];
rz(-0.94347763) q[3];
sx q[3];
rz(0.77997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778377) q[0];
sx q[0];
rz(-1.5564593) q[0];
sx q[0];
rz(-2.7177287) q[0];
rz(-2.9546812) q[1];
sx q[1];
rz(-0.82004768) q[1];
sx q[1];
rz(-1.6835469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20987147) q[0];
sx q[0];
rz(-0.11882028) q[0];
sx q[0];
rz(-3.0221536) q[0];
rz(-pi) q[1];
rz(-2.0844314) q[2];
sx q[2];
rz(-1.8087401) q[2];
sx q[2];
rz(0.4713716) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9002237) q[1];
sx q[1];
rz(-1.4321064) q[1];
sx q[1];
rz(-1.6226416) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23641674) q[3];
sx q[3];
rz(-2.7763753) q[3];
sx q[3];
rz(-1.055298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75108782) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(-2.0474153) q[2];
rz(1.4607726) q[3];
sx q[3];
rz(-1.3760309) q[3];
sx q[3];
rz(-1.338753) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8193034) q[0];
sx q[0];
rz(-1.3811454) q[0];
sx q[0];
rz(1.5018916) q[0];
rz(0.27885258) q[1];
sx q[1];
rz(-2.1809705) q[1];
sx q[1];
rz(1.0241114) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8968643) q[0];
sx q[0];
rz(-0.49806589) q[0];
sx q[0];
rz(-0.19356541) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7488519) q[2];
sx q[2];
rz(-2.8626056) q[2];
sx q[2];
rz(-1.6131608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1361724) q[1];
sx q[1];
rz(-1.3498303) q[1];
sx q[1];
rz(-1.7762089) q[1];
x q[2];
rz(-0.34124438) q[3];
sx q[3];
rz(-2.2856718) q[3];
sx q[3];
rz(1.224186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2037105) q[2];
sx q[2];
rz(-1.8546162) q[2];
sx q[2];
rz(2.6046275) q[2];
rz(1.9133866) q[3];
sx q[3];
rz(-2.1824586) q[3];
sx q[3];
rz(-2.8502407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1163597) q[0];
sx q[0];
rz(-1.3675084) q[0];
sx q[0];
rz(1.0687923) q[0];
rz(-0.7069201) q[1];
sx q[1];
rz(-1.128935) q[1];
sx q[1];
rz(3.1405906) q[1];
rz(2.7189485) q[2];
sx q[2];
rz(-0.44036897) q[2];
sx q[2];
rz(-1.8457495) q[2];
rz(-0.14866004) q[3];
sx q[3];
rz(-1.3547263) q[3];
sx q[3];
rz(2.3480036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
