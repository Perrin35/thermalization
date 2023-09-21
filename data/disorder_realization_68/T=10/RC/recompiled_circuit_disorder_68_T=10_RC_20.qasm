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
rz(-0.49180254) q[0];
sx q[0];
rz(-2.9536182) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(4.6586577) q[1];
sx q[1];
rz(12.933856) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7029019) q[0];
sx q[0];
rz(-0.54649788) q[0];
sx q[0];
rz(1.7706857) q[0];
rz(-pi) q[1];
rz(1.3221402) q[2];
sx q[2];
rz(-0.50422943) q[2];
sx q[2];
rz(2.3766975) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5012813) q[1];
sx q[1];
rz(-0.2959364) q[1];
sx q[1];
rz(-0.96247767) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7484619) q[3];
sx q[3];
rz(-1.9068309) q[3];
sx q[3];
rz(-2.6538268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.964103) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-2.5906079) q[2];
rz(1.8356813) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(-1.8252385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6630163) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(-2.7117803) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(-2.205251) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4124356) q[0];
sx q[0];
rz(-1.4834187) q[0];
sx q[0];
rz(-2.9108414) q[0];
rz(-2.721644) q[2];
sx q[2];
rz(-2.311085) q[2];
sx q[2];
rz(0.74479693) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71948235) q[1];
sx q[1];
rz(-1.6750095) q[1];
sx q[1];
rz(-0.69394333) q[1];
rz(1.6132658) q[3];
sx q[3];
rz(-2.6385033) q[3];
sx q[3];
rz(-0.79723061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77461809) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(-0.42638391) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(-3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8957829) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(-0.93908969) q[0];
rz(0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(0.59392196) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1547326) q[0];
sx q[0];
rz(-1.2504471) q[0];
sx q[0];
rz(-0.31517584) q[0];
x q[1];
rz(-1.6410286) q[2];
sx q[2];
rz(-0.68968455) q[2];
sx q[2];
rz(2.1948512) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18332874) q[1];
sx q[1];
rz(-1.3576344) q[1];
sx q[1];
rz(2.0002685) q[1];
x q[2];
rz(-1.2265615) q[3];
sx q[3];
rz(-1.1262745) q[3];
sx q[3];
rz(-1.759474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64017355) q[2];
sx q[2];
rz(-0.72945014) q[2];
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
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664292) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(-2.638812) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(-0.75685135) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3319791) q[0];
sx q[0];
rz(-2.6410714) q[0];
sx q[0];
rz(-2.39397) q[0];
x q[1];
rz(1.1050622) q[2];
sx q[2];
rz(-1.5823936) q[2];
sx q[2];
rz(0.84601814) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6824324) q[1];
sx q[1];
rz(-1.2130514) q[1];
sx q[1];
rz(-0.95174241) q[1];
x q[2];
rz(2.5960856) q[3];
sx q[3];
rz(-2.3295998) q[3];
sx q[3];
rz(-3.0363887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.42671529) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(1.4871917) q[2];
rz(-0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(-2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-0.75772444) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(1.0505189) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3462853) q[0];
sx q[0];
rz(-1.8827569) q[0];
sx q[0];
rz(0.26766582) q[0];
rz(-0.97425766) q[2];
sx q[2];
rz(-1.9447118) q[2];
sx q[2];
rz(0.24335441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95549059) q[1];
sx q[1];
rz(-2.0932066) q[1];
sx q[1];
rz(2.2239457) q[1];
rz(-pi) q[2];
rz(-0.80661185) q[3];
sx q[3];
rz(-1.9730554) q[3];
sx q[3];
rz(2.3276687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.918255) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(2.5081432) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(-1.6075915) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(-1.4621428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7397241) q[0];
sx q[0];
rz(-1.7812294) q[0];
sx q[0];
rz(1.4896638) q[0];
x q[1];
rz(2.9003733) q[2];
sx q[2];
rz(-0.93572817) q[2];
sx q[2];
rz(-1.5347753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0837005) q[1];
sx q[1];
rz(-1.18827) q[1];
sx q[1];
rz(1.7544569) q[1];
rz(-1.2268279) q[3];
sx q[3];
rz(-1.2729537) q[3];
sx q[3];
rz(-2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2016466) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(-2.3366826) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686304) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(-2.1169128) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(1.0120846) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1597848) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(2.3083789) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7527533) q[2];
sx q[2];
rz(-1.6992339) q[2];
sx q[2];
rz(1.5156137) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95841366) q[1];
sx q[1];
rz(-1.603754) q[1];
sx q[1];
rz(2.5977913) q[1];
rz(-pi) q[2];
rz(0.55862553) q[3];
sx q[3];
rz(-1.6405676) q[3];
sx q[3];
rz(0.39486265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-2.7116595) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5291418) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(2.9274143) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(2.8578551) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3389694) q[0];
sx q[0];
rz(-1.6049275) q[0];
sx q[0];
rz(-2.8404833) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9494638) q[2];
sx q[2];
rz(-1.8953952) q[2];
sx q[2];
rz(-1.0629551) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6730496) q[1];
sx q[1];
rz(-1.0769516) q[1];
sx q[1];
rz(-1.7935351) q[1];
rz(-pi) q[2];
x q[2];
rz(2.22105) q[3];
sx q[3];
rz(-2.1170108) q[3];
sx q[3];
rz(1.2342681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6909137) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(-1.4452176) q[2];
rz(-1.5444267) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504869) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(-3.0723363) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(-1.5725296) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5208961) q[0];
sx q[0];
rz(-2.3224152) q[0];
sx q[0];
rz(2.2197414) q[0];
rz(1.6530767) q[2];
sx q[2];
rz(-1.2459718) q[2];
sx q[2];
rz(-2.5282853) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3410586) q[1];
sx q[1];
rz(-1.8039928) q[1];
sx q[1];
rz(-1.484616) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92215718) q[3];
sx q[3];
rz(-1.7382009) q[3];
sx q[3];
rz(2.5006014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1853603) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(3.0017079) q[2];
rz(0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(1.2311252) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(2.8737601) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56775996) q[0];
sx q[0];
rz(-2.0127675) q[0];
sx q[0];
rz(-0.7451591) q[0];
x q[1];
rz(1.7540625) q[2];
sx q[2];
rz(-1.4472618) q[2];
sx q[2];
rz(-1.087041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1416727) q[1];
sx q[1];
rz(-1.281764) q[1];
sx q[1];
rz(0.33860597) q[1];
rz(-2.855741) q[3];
sx q[3];
rz(-2.0945858) q[3];
sx q[3];
rz(-1.3956192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.81007593) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.4769185) q[2];
rz(-2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(0.5464856) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(-2.3616882) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(-0.65200381) q[2];
sx q[2];
rz(-2.3163788) q[2];
sx q[2];
rz(-0.083995081) q[2];
rz(-0.6775425) q[3];
sx q[3];
rz(-1.7384221) q[3];
sx q[3];
rz(1.3330028) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];