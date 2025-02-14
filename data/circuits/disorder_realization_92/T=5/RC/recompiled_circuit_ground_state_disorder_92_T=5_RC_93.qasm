OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.3761223) q[0];
sx q[0];
rz(2.0630615) q[0];
sx q[0];
rz(10.800092) q[0];
rz(-2.7219661) q[1];
sx q[1];
rz(-1.3803177) q[1];
sx q[1];
rz(-0.79545155) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8859098) q[0];
sx q[0];
rz(-1.1982746) q[0];
sx q[0];
rz(1.4857948) q[0];
x q[1];
rz(1.4515204) q[2];
sx q[2];
rz(-1.1732427) q[2];
sx q[2];
rz(1.9089323) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1297721) q[1];
sx q[1];
rz(-1.3737965) q[1];
sx q[1];
rz(-2.4596766) q[1];
rz(1.8849089) q[3];
sx q[3];
rz(-2.8693349) q[3];
sx q[3];
rz(1.5262878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85806435) q[2];
sx q[2];
rz(-0.29164803) q[2];
sx q[2];
rz(-2.7533599) q[2];
rz(-0.20902108) q[3];
sx q[3];
rz(-1.6671323) q[3];
sx q[3];
rz(2.2504263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46880576) q[0];
sx q[0];
rz(-2.6728632) q[0];
sx q[0];
rz(2.3720473) q[0];
rz(-1.0308713) q[1];
sx q[1];
rz(-2.2831235) q[1];
sx q[1];
rz(1.7147725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1772004) q[0];
sx q[0];
rz(-1.3528723) q[0];
sx q[0];
rz(0.30857463) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.513711) q[2];
sx q[2];
rz(-2.2453614) q[2];
sx q[2];
rz(0.13007535) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64411341) q[1];
sx q[1];
rz(-2.1783324) q[1];
sx q[1];
rz(1.8827022) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2852422) q[3];
sx q[3];
rz(-1.5305686) q[3];
sx q[3];
rz(2.2742607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1398805) q[2];
sx q[2];
rz(-2.6931245) q[2];
sx q[2];
rz(1.2343538) q[2];
rz(2.6185696) q[3];
sx q[3];
rz(-1.1250027) q[3];
sx q[3];
rz(2.4856868) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6958385) q[0];
sx q[0];
rz(-0.0025175968) q[0];
sx q[0];
rz(-0.54760951) q[0];
rz(0.86946636) q[1];
sx q[1];
rz(-2.165386) q[1];
sx q[1];
rz(-1.4617823) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6637409) q[0];
sx q[0];
rz(-1.9123075) q[0];
sx q[0];
rz(2.6966901) q[0];
rz(-pi) q[1];
rz(1.3140979) q[2];
sx q[2];
rz(-0.82523286) q[2];
sx q[2];
rz(-1.0302533) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28785767) q[1];
sx q[1];
rz(-2.6035941) q[1];
sx q[1];
rz(-0.96189672) q[1];
rz(2.2814155) q[3];
sx q[3];
rz(-1.0499991) q[3];
sx q[3];
rz(1.8627244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.72472858) q[2];
sx q[2];
rz(-2.8432196) q[2];
sx q[2];
rz(0.63428632) q[2];
rz(2.0554845) q[3];
sx q[3];
rz(-0.48091286) q[3];
sx q[3];
rz(-1.1336795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61741614) q[0];
sx q[0];
rz(-2.9845181) q[0];
sx q[0];
rz(-1.5217391) q[0];
rz(2.1583083) q[1];
sx q[1];
rz(-2.0359928) q[1];
sx q[1];
rz(-0.081667893) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3043609) q[0];
sx q[0];
rz(-1.0153234) q[0];
sx q[0];
rz(0.13120478) q[0];
rz(1.2377023) q[2];
sx q[2];
rz(-1.5577847) q[2];
sx q[2];
rz(-1.5481373) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5809243) q[1];
sx q[1];
rz(-0.67519516) q[1];
sx q[1];
rz(-0.84431727) q[1];
rz(1.6713167) q[3];
sx q[3];
rz(-0.49197061) q[3];
sx q[3];
rz(0.63382404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4982768) q[2];
sx q[2];
rz(-1.4772819) q[2];
sx q[2];
rz(2.2225883) q[2];
rz(-0.32215858) q[3];
sx q[3];
rz(-1.5124269) q[3];
sx q[3];
rz(-2.7951516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084694013) q[0];
sx q[0];
rz(-2.0094805) q[0];
sx q[0];
rz(-1.3380916) q[0];
rz(0.40254205) q[1];
sx q[1];
rz(-1.891338) q[1];
sx q[1];
rz(-1.5571627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8020222) q[0];
sx q[0];
rz(-2.5568107) q[0];
sx q[0];
rz(1.0735548) q[0];
rz(-1.8025777) q[2];
sx q[2];
rz(-2.351481) q[2];
sx q[2];
rz(-2.6274644) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8338523) q[1];
sx q[1];
rz(-0.55451143) q[1];
sx q[1];
rz(-0.28077287) q[1];
rz(2.8731737) q[3];
sx q[3];
rz(-1.1158593) q[3];
sx q[3];
rz(1.1283629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1160139) q[2];
sx q[2];
rz(-1.8639001) q[2];
sx q[2];
rz(-0.73520994) q[2];
rz(2.0871346) q[3];
sx q[3];
rz(-1.0337044) q[3];
sx q[3];
rz(1.2247156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16578199) q[0];
sx q[0];
rz(-2.7175792) q[0];
sx q[0];
rz(-1.2531511) q[0];
rz(-1.5660628) q[1];
sx q[1];
rz(-1.1352481) q[1];
sx q[1];
rz(-2.6083686) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76372242) q[0];
sx q[0];
rz(-1.8061751) q[0];
sx q[0];
rz(-0.81416582) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27006939) q[2];
sx q[2];
rz(-0.52723072) q[2];
sx q[2];
rz(-0.6747077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27665177) q[1];
sx q[1];
rz(-1.3485695) q[1];
sx q[1];
rz(-1.3781616) q[1];
x q[2];
rz(0.56508692) q[3];
sx q[3];
rz(-2.0425519) q[3];
sx q[3];
rz(0.12607546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.48625654) q[2];
sx q[2];
rz(-0.37898263) q[2];
sx q[2];
rz(0.8302702) q[2];
rz(-1.8298979) q[3];
sx q[3];
rz(-2.0528841) q[3];
sx q[3];
rz(-0.24454288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.526392) q[0];
sx q[0];
rz(-1.4854234) q[0];
sx q[0];
rz(-2.3080589) q[0];
rz(2.1444881) q[1];
sx q[1];
rz(-1.763696) q[1];
sx q[1];
rz(-1.5514143) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1520321) q[0];
sx q[0];
rz(-0.44267926) q[0];
sx q[0];
rz(0.17501207) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0710414) q[2];
sx q[2];
rz(-0.65325551) q[2];
sx q[2];
rz(-0.61966255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.60234705) q[1];
sx q[1];
rz(-1.6825469) q[1];
sx q[1];
rz(1.4503327) q[1];
rz(-pi) q[2];
rz(2.6069712) q[3];
sx q[3];
rz(-1.3546913) q[3];
sx q[3];
rz(-0.099319746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3414063) q[2];
sx q[2];
rz(-2.1087346) q[2];
sx q[2];
rz(1.3979073) q[2];
rz(0.37311113) q[3];
sx q[3];
rz(-2.2171376) q[3];
sx q[3];
rz(1.0188518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7031192) q[0];
sx q[0];
rz(-1.4716453) q[0];
sx q[0];
rz(2.8224509) q[0];
rz(1.1646264) q[1];
sx q[1];
rz(-1.9705557) q[1];
sx q[1];
rz(-3.107531) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1036164) q[0];
sx q[0];
rz(-1.8744133) q[0];
sx q[0];
rz(-2.3453317) q[0];
rz(-pi) q[1];
rz(-2.3715842) q[2];
sx q[2];
rz(-1.7988552) q[2];
sx q[2];
rz(2.8213825) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0314591) q[1];
sx q[1];
rz(-1.5601399) q[1];
sx q[1];
rz(1.5489122) q[1];
rz(1.930792) q[3];
sx q[3];
rz(-1.0673095) q[3];
sx q[3];
rz(0.22758037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.065319149) q[2];
sx q[2];
rz(-1.3955782) q[2];
sx q[2];
rz(0.070153959) q[2];
rz(-1.8326727) q[3];
sx q[3];
rz(-2.2891243) q[3];
sx q[3];
rz(0.13516983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081130505) q[0];
sx q[0];
rz(-1.9996996) q[0];
sx q[0];
rz(0.4749701) q[0];
rz(1.9604856) q[1];
sx q[1];
rz(-0.89473748) q[1];
sx q[1];
rz(2.2244577) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346996) q[0];
sx q[0];
rz(-2.1297161) q[0];
sx q[0];
rz(-1.8447645) q[0];
rz(1.195329) q[2];
sx q[2];
rz(-2.4175543) q[2];
sx q[2];
rz(2.6065116) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5702906) q[1];
sx q[1];
rz(-1.362974) q[1];
sx q[1];
rz(0.083870693) q[1];
rz(-pi) q[2];
rz(1.9539509) q[3];
sx q[3];
rz(-2.1360548) q[3];
sx q[3];
rz(1.7704028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8721623) q[2];
sx q[2];
rz(-1.3698801) q[2];
sx q[2];
rz(-1.5577462) q[2];
rz(2.9861084) q[3];
sx q[3];
rz(-2.1126066) q[3];
sx q[3];
rz(2.1492929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48142049) q[0];
sx q[0];
rz(-1.6654797) q[0];
sx q[0];
rz(0.36556622) q[0];
rz(1.9407326) q[1];
sx q[1];
rz(-1.6114085) q[1];
sx q[1];
rz(-1.0467122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4037673) q[0];
sx q[0];
rz(-1.0407454) q[0];
sx q[0];
rz(0.60204317) q[0];
rz(0.44136845) q[2];
sx q[2];
rz(-1.0689931) q[2];
sx q[2];
rz(-0.29870118) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9871708) q[1];
sx q[1];
rz(-0.94564309) q[1];
sx q[1];
rz(-0.67657569) q[1];
rz(-pi) q[2];
rz(0.58308122) q[3];
sx q[3];
rz(-2.0399391) q[3];
sx q[3];
rz(-0.52419477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6863579) q[2];
sx q[2];
rz(-1.1150259) q[2];
sx q[2];
rz(-2.4380747) q[2];
rz(2.5840058) q[3];
sx q[3];
rz(-2.360354) q[3];
sx q[3];
rz(-2.668146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8483509) q[0];
sx q[0];
rz(-2.1108755) q[0];
sx q[0];
rz(-2.1847771) q[0];
rz(2.1059857) q[1];
sx q[1];
rz(-2.5685812) q[1];
sx q[1];
rz(-2.0249637) q[1];
rz(0.63450068) q[2];
sx q[2];
rz(-2.366323) q[2];
sx q[2];
rz(0.97699588) q[2];
rz(-1.5286636) q[3];
sx q[3];
rz(-1.7880472) q[3];
sx q[3];
rz(0.40732297) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
