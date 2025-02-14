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
rz(0.19658495) q[0];
sx q[0];
rz(-2.6994446) q[0];
sx q[0];
rz(2.2801939) q[0];
rz(-1.6317033) q[1];
sx q[1];
rz(-1.3246526) q[1];
sx q[1];
rz(-3.1060001) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7073878) q[0];
sx q[0];
rz(-1.2363096) q[0];
sx q[0];
rz(-1.3637961) q[0];
x q[1];
rz(-1.6794559) q[2];
sx q[2];
rz(-1.265996) q[2];
sx q[2];
rz(-2.5632312) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10251789) q[1];
sx q[1];
rz(-0.84474428) q[1];
sx q[1];
rz(1.6186569) q[1];
rz(-pi) q[2];
rz(2.2507812) q[3];
sx q[3];
rz(-1.5847795) q[3];
sx q[3];
rz(0.76410642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5888136) q[2];
sx q[2];
rz(-1.3014883) q[2];
sx q[2];
rz(2.0471052) q[2];
rz(1.5158481) q[3];
sx q[3];
rz(-0.87259126) q[3];
sx q[3];
rz(2.998013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0933519) q[0];
sx q[0];
rz(-0.99026647) q[0];
sx q[0];
rz(2.7675203) q[0];
rz(0.85809842) q[1];
sx q[1];
rz(-0.94081196) q[1];
sx q[1];
rz(-2.0657952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25903048) q[0];
sx q[0];
rz(-0.1452862) q[0];
sx q[0];
rz(2.3502716) q[0];
x q[1];
rz(1.3316989) q[2];
sx q[2];
rz(-1.0361203) q[2];
sx q[2];
rz(2.3102674) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.058671) q[1];
sx q[1];
rz(-1.4545014) q[1];
sx q[1];
rz(2.7398159) q[1];
x q[2];
rz(0.91173633) q[3];
sx q[3];
rz(-1.0056595) q[3];
sx q[3];
rz(-0.46253461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51521236) q[2];
sx q[2];
rz(-0.43312803) q[2];
sx q[2];
rz(-2.058378) q[2];
rz(1.8488688) q[3];
sx q[3];
rz(-1.1336528) q[3];
sx q[3];
rz(2.2129272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.885963) q[0];
sx q[0];
rz(-0.77115458) q[0];
sx q[0];
rz(-2.6166925) q[0];
rz(1.6820172) q[1];
sx q[1];
rz(-2.4039905) q[1];
sx q[1];
rz(2.9958013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12851579) q[0];
sx q[0];
rz(-1.4798963) q[0];
sx q[0];
rz(-0.10581067) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1728233) q[2];
sx q[2];
rz(-1.1175691) q[2];
sx q[2];
rz(-2.9211958) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3891366) q[1];
sx q[1];
rz(-2.698878) q[1];
sx q[1];
rz(2.929844) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1946241) q[3];
sx q[3];
rz(-2.3071529) q[3];
sx q[3];
rz(-0.26716993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8968481) q[2];
sx q[2];
rz(-1.1620099) q[2];
sx q[2];
rz(-2.9761918) q[2];
rz(-2.2798955) q[3];
sx q[3];
rz(-2.4498037) q[3];
sx q[3];
rz(-2.944788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8621314) q[0];
sx q[0];
rz(-2.617351) q[0];
sx q[0];
rz(0.72823802) q[0];
rz(-0.32314745) q[1];
sx q[1];
rz(-2.1073982) q[1];
sx q[1];
rz(1.039215) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0511) q[0];
sx q[0];
rz(-1.1742102) q[0];
sx q[0];
rz(-2.0329679) q[0];
x q[1];
rz(2.8604554) q[2];
sx q[2];
rz(-2.1875811) q[2];
sx q[2];
rz(-2.5532818) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2715214) q[1];
sx q[1];
rz(-2.0539375) q[1];
sx q[1];
rz(1.6752233) q[1];
rz(2.3023241) q[3];
sx q[3];
rz(-0.89863649) q[3];
sx q[3];
rz(1.9306926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6497961) q[2];
sx q[2];
rz(-2.0911262) q[2];
sx q[2];
rz(0.51898471) q[2];
rz(2.0670048) q[3];
sx q[3];
rz(-1.855775) q[3];
sx q[3];
rz(0.49526596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0786667) q[0];
sx q[0];
rz(-2.2195897) q[0];
sx q[0];
rz(-0.17157383) q[0];
rz(-1.0942787) q[1];
sx q[1];
rz(-1.3010052) q[1];
sx q[1];
rz(2.9069854) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0820935) q[0];
sx q[0];
rz(-1.4429868) q[0];
sx q[0];
rz(-1.2907842) q[0];
rz(-pi) q[1];
rz(-1.6233088) q[2];
sx q[2];
rz(-1.0641452) q[2];
sx q[2];
rz(1.5280629) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68069862) q[1];
sx q[1];
rz(-1.0061703) q[1];
sx q[1];
rz(2.2959397) q[1];
x q[2];
rz(0.42886491) q[3];
sx q[3];
rz(-2.3842616) q[3];
sx q[3];
rz(-1.1348108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2285063) q[2];
sx q[2];
rz(-1.191491) q[2];
sx q[2];
rz(-1.7388434) q[2];
rz(-0.62471041) q[3];
sx q[3];
rz(-0.96840817) q[3];
sx q[3];
rz(-1.3431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45873555) q[0];
sx q[0];
rz(-3.1338437) q[0];
sx q[0];
rz(-0.32753456) q[0];
rz(-0.028118357) q[1];
sx q[1];
rz(-1.997812) q[1];
sx q[1];
rz(-0.99172529) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1338746) q[0];
sx q[0];
rz(-1.0622866) q[0];
sx q[0];
rz(2.4009812) q[0];
x q[1];
rz(1.6794793) q[2];
sx q[2];
rz(-1.7955488) q[2];
sx q[2];
rz(-1.0959742) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8649329) q[1];
sx q[1];
rz(-2.3845379) q[1];
sx q[1];
rz(-1.6318342) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5434274) q[3];
sx q[3];
rz(-1.1002161) q[3];
sx q[3];
rz(-1.9675627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.907054) q[2];
sx q[2];
rz(-1.3776642) q[2];
sx q[2];
rz(-1.482359) q[2];
rz(2.0756857) q[3];
sx q[3];
rz(-1.1803455) q[3];
sx q[3];
rz(-2.0445686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3858353) q[0];
sx q[0];
rz(-3.0297854) q[0];
sx q[0];
rz(-0.29543153) q[0];
rz(2.9040728) q[1];
sx q[1];
rz(-2.2078881) q[1];
sx q[1];
rz(2.3942153) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4656671) q[0];
sx q[0];
rz(-0.79361483) q[0];
sx q[0];
rz(2.3966952) q[0];
rz(2.105999) q[2];
sx q[2];
rz(-1.0135302) q[2];
sx q[2];
rz(-2.0095428) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.566639) q[1];
sx q[1];
rz(-2.5051834) q[1];
sx q[1];
rz(-2.2242935) q[1];
rz(-pi) q[2];
rz(-2.5980972) q[3];
sx q[3];
rz(-1.8955909) q[3];
sx q[3];
rz(-0.80431688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6682917) q[2];
sx q[2];
rz(-1.95582) q[2];
sx q[2];
rz(-0.2549003) q[2];
rz(-2.9292246) q[3];
sx q[3];
rz(-2.2771211) q[3];
sx q[3];
rz(0.00096850639) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2905228) q[0];
sx q[0];
rz(-1.9177328) q[0];
sx q[0];
rz(3.0176924) q[0];
rz(-2.2663785) q[1];
sx q[1];
rz(-0.83299914) q[1];
sx q[1];
rz(2.5828054) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4717455) q[0];
sx q[0];
rz(-0.59006834) q[0];
sx q[0];
rz(-1.1099932) q[0];
x q[1];
rz(-0.28382878) q[2];
sx q[2];
rz(-2.0528194) q[2];
sx q[2];
rz(1.986077) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81441036) q[1];
sx q[1];
rz(-0.60403999) q[1];
sx q[1];
rz(-2.6534897) q[1];
x q[2];
rz(1.0304006) q[3];
sx q[3];
rz(-2.4321788) q[3];
sx q[3];
rz(-1.7075001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.045804068) q[2];
sx q[2];
rz(-0.49171058) q[2];
sx q[2];
rz(2.4208505) q[2];
rz(-0.36302429) q[3];
sx q[3];
rz(-1.2237153) q[3];
sx q[3];
rz(-1.5117517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95661288) q[0];
sx q[0];
rz(-2.9473801) q[0];
sx q[0];
rz(-0.089056253) q[0];
rz(-0.36674276) q[1];
sx q[1];
rz(-0.9042424) q[1];
sx q[1];
rz(1.6820224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7686667) q[0];
sx q[0];
rz(-2.4690876) q[0];
sx q[0];
rz(-1.1993394) q[0];
x q[1];
rz(-0.681293) q[2];
sx q[2];
rz(-2.3325787) q[2];
sx q[2];
rz(1.9495277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68946099) q[1];
sx q[1];
rz(-2.2575827) q[1];
sx q[1];
rz(0.48312123) q[1];
rz(-pi) q[2];
rz(-2.5427548) q[3];
sx q[3];
rz(-1.2019079) q[3];
sx q[3];
rz(1.2303331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2064994) q[2];
sx q[2];
rz(-1.7743899) q[2];
sx q[2];
rz(1.8915141) q[2];
rz(1.2262723) q[3];
sx q[3];
rz(-0.63392249) q[3];
sx q[3];
rz(2.5026076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117689) q[0];
sx q[0];
rz(-2.0084232) q[0];
sx q[0];
rz(2.0015707) q[0];
rz(-0.58311588) q[1];
sx q[1];
rz(-1.6551599) q[1];
sx q[1];
rz(-2.4737632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4224189) q[0];
sx q[0];
rz(-1.2454709) q[0];
sx q[0];
rz(-2.3803985) q[0];
x q[1];
rz(2.2726353) q[2];
sx q[2];
rz(-1.387082) q[2];
sx q[2];
rz(-0.76317235) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9676303) q[1];
sx q[1];
rz(-0.47904992) q[1];
sx q[1];
rz(1.1381989) q[1];
x q[2];
rz(-2.3488914) q[3];
sx q[3];
rz(-1.2374733) q[3];
sx q[3];
rz(-1.9938716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2976611) q[2];
sx q[2];
rz(-0.49489489) q[2];
sx q[2];
rz(-2.3893791) q[2];
rz(3.0609868) q[3];
sx q[3];
rz(-1.3496496) q[3];
sx q[3];
rz(2.5940671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053631) q[0];
sx q[0];
rz(-1.6148051) q[0];
sx q[0];
rz(-0.057407277) q[0];
rz(-2.2930131) q[1];
sx q[1];
rz(-0.72863693) q[1];
sx q[1];
rz(-1.4562664) q[1];
rz(-2.3586629) q[2];
sx q[2];
rz(-0.83275262) q[2];
sx q[2];
rz(0.61093753) q[2];
rz(-0.46075321) q[3];
sx q[3];
rz(-1.2450153) q[3];
sx q[3];
rz(1.0871441) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
