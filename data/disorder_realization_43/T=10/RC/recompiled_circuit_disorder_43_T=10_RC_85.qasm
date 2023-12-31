OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(1.8068846) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0358409) q[0];
sx q[0];
rz(-0.63397206) q[0];
sx q[0];
rz(-2.7266154) q[0];
rz(-pi) q[1];
rz(-2.072233) q[2];
sx q[2];
rz(-0.62383365) q[2];
sx q[2];
rz(1.2944348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18260278) q[1];
sx q[1];
rz(-2.5187153) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(1.9777771) q[3];
sx q[3];
rz(-2.3325936) q[3];
sx q[3];
rz(-2.80796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(-1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48288229) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(-3.1352502) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77942383) q[0];
sx q[0];
rz(-0.2340901) q[0];
sx q[0];
rz(0.97658821) q[0];
rz(2.2750521) q[2];
sx q[2];
rz(-0.73359493) q[2];
sx q[2];
rz(1.7689592) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.43286846) q[1];
sx q[1];
rz(-1.2607062) q[1];
sx q[1];
rz(-1.5061782) q[1];
rz(-2.7553495) q[3];
sx q[3];
rz(-2.5442903) q[3];
sx q[3];
rz(-0.18988767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(0.0058962065) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5851615) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(-2.7745461) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(-0.095741622) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75671065) q[0];
sx q[0];
rz(-1.7236992) q[0];
sx q[0];
rz(1.3923313) q[0];
rz(-0.7272561) q[2];
sx q[2];
rz(-0.63823344) q[2];
sx q[2];
rz(-2.9122695) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9853471) q[1];
sx q[1];
rz(-0.74838446) q[1];
sx q[1];
rz(-2.2845539) q[1];
x q[2];
rz(-1.5355574) q[3];
sx q[3];
rz(-2.1217151) q[3];
sx q[3];
rz(-0.86722022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0720955) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(-0.21162027) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(0.37477469) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(1.0338763) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046712) q[0];
sx q[0];
rz(-1.0949507) q[0];
sx q[0];
rz(-0.67142077) q[0];
rz(2.3029762) q[2];
sx q[2];
rz(-0.71574434) q[2];
sx q[2];
rz(-1.0775623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4775866) q[1];
sx q[1];
rz(-2.9395182) q[1];
sx q[1];
rz(1.2408153) q[1];
rz(-pi) q[2];
rz(-0.86446188) q[3];
sx q[3];
rz(-0.41941386) q[3];
sx q[3];
rz(2.1692587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.7948077) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(2.72686) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(2.5674852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60526472) q[0];
sx q[0];
rz(-0.61075532) q[0];
sx q[0];
rz(1.3461793) q[0];
rz(-pi) q[1];
rz(-2.4262869) q[2];
sx q[2];
rz(-1.5537405) q[2];
sx q[2];
rz(-2.7280083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.638006) q[1];
sx q[1];
rz(-1.9203016) q[1];
sx q[1];
rz(-1.867678) q[1];
rz(-2.2361034) q[3];
sx q[3];
rz(-1.724616) q[3];
sx q[3];
rz(0.32115768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85270143) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.627702) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742652) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(0.69818991) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(2.4093157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9376611) q[0];
sx q[0];
rz(-1.7294149) q[0];
sx q[0];
rz(-2.926814) q[0];
x q[1];
rz(2.4562624) q[2];
sx q[2];
rz(-2.2517859) q[2];
sx q[2];
rz(-0.53390098) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60019833) q[1];
sx q[1];
rz(-1.3458034) q[1];
sx q[1];
rz(0.66585983) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3527649) q[3];
sx q[3];
rz(-2.1302967) q[3];
sx q[3];
rz(-0.41527173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7531062) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.1674315) q[3];
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
rz(0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(0.15643315) q[0];
rz(-0.45174831) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(2.3715473) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5891588) q[0];
sx q[0];
rz(-1.3981515) q[0];
sx q[0];
rz(0.60683672) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1427878) q[2];
sx q[2];
rz(-0.6711798) q[2];
sx q[2];
rz(-2.3064248) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60939497) q[1];
sx q[1];
rz(-2.1708793) q[1];
sx q[1];
rz(2.1983912) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5823334) q[3];
sx q[3];
rz(-2.5659605) q[3];
sx q[3];
rz(1.5534793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4975171) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(1.0151803) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-2.6530478) q[1];
sx q[1];
rz(0.36639211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05224932) q[0];
sx q[0];
rz(-1.0567259) q[0];
sx q[0];
rz(0.71787562) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5744152) q[2];
sx q[2];
rz(-0.71225538) q[2];
sx q[2];
rz(0.30356193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39216343) q[1];
sx q[1];
rz(-2.2656419) q[1];
sx q[1];
rz(1.7872582) q[1];
rz(2.1226235) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(-3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(0.29385847) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(-2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(0.88395399) q[0];
rz(2.4813095) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(-2.3502137) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34245472) q[0];
sx q[0];
rz(-0.90478169) q[0];
sx q[0];
rz(2.7333583) q[0];
x q[1];
rz(-1.0857401) q[2];
sx q[2];
rz(-1.9888478) q[2];
sx q[2];
rz(-2.5784091) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8225122) q[1];
sx q[1];
rz(-1.1573536) q[1];
sx q[1];
rz(2.6507069) q[1];
rz(0.16718849) q[3];
sx q[3];
rz(-1.0382004) q[3];
sx q[3];
rz(-2.5933468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(2.1203314) q[2];
rz(2.7630473) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1148949) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(2.5337906) q[0];
rz(-0.23884493) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(1.3806608) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138807) q[0];
sx q[0];
rz(-0.6427592) q[0];
sx q[0];
rz(1.3520665) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5153377) q[2];
sx q[2];
rz(-0.43606731) q[2];
sx q[2];
rz(0.01576327) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7541411) q[1];
sx q[1];
rz(-1.8879461) q[1];
sx q[1];
rz(-1.1412568) q[1];
rz(-pi) q[2];
rz(-0.72762604) q[3];
sx q[3];
rz(-2.8907223) q[3];
sx q[3];
rz(0.48654702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-0.33622462) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(2.5972988) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(0.33967321) q[2];
sx q[2];
rz(-0.85396955) q[2];
sx q[2];
rz(0.11243482) q[2];
rz(1.929677) q[3];
sx q[3];
rz(-2.5419895) q[3];
sx q[3];
rz(0.06751577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
