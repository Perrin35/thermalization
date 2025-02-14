OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(2.03696) q[0];
rz(-1.8072577) q[1];
sx q[1];
rz(5.5601064) q[1];
sx q[1];
rz(13.830166) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4975464) q[0];
sx q[0];
rz(-0.9263557) q[0];
sx q[0];
rz(-0.30797663) q[0];
x q[1];
rz(-0.99743263) q[2];
sx q[2];
rz(-1.9321529) q[2];
sx q[2];
rz(-0.45623764) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3448779) q[1];
sx q[1];
rz(-1.895322) q[1];
sx q[1];
rz(3.0561563) q[1];
rz(2.0094517) q[3];
sx q[3];
rz(-1.3490455) q[3];
sx q[3];
rz(-0.36166276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9870712) q[2];
sx q[2];
rz(-1.0769341) q[2];
sx q[2];
rz(-1.7729574) q[2];
rz(-2.9030419) q[3];
sx q[3];
rz(-0.41540256) q[3];
sx q[3];
rz(-3.0273048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7411165) q[0];
sx q[0];
rz(-1.1392765) q[0];
sx q[0];
rz(1.9912632) q[0];
rz(0.087609619) q[1];
sx q[1];
rz(-1.3299512) q[1];
sx q[1];
rz(-1.5706583) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.621839) q[0];
sx q[0];
rz(-1.6770419) q[0];
sx q[0];
rz(0.48522093) q[0];
rz(2.9628062) q[2];
sx q[2];
rz(-1.6331722) q[2];
sx q[2];
rz(-1.5449926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1521235) q[1];
sx q[1];
rz(-1.996576) q[1];
sx q[1];
rz(-0.87321891) q[1];
rz(-pi) q[2];
rz(0.79227607) q[3];
sx q[3];
rz(-1.918692) q[3];
sx q[3];
rz(-2.874305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7972083) q[2];
sx q[2];
rz(-0.71343652) q[2];
sx q[2];
rz(0.1864645) q[2];
rz(-0.79948419) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(-2.1605087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8568273) q[0];
sx q[0];
rz(-1.3815877) q[0];
sx q[0];
rz(3.0322266) q[0];
rz(0.60375396) q[1];
sx q[1];
rz(-1.8021288) q[1];
sx q[1];
rz(1.0341136) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4277735) q[0];
sx q[0];
rz(-1.7805011) q[0];
sx q[0];
rz(-0.13596491) q[0];
x q[1];
rz(-1.6414406) q[2];
sx q[2];
rz(-1.4817258) q[2];
sx q[2];
rz(-1.3180817) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2769985) q[1];
sx q[1];
rz(-1.8499476) q[1];
sx q[1];
rz(1.6807081) q[1];
rz(-1.7121186) q[3];
sx q[3];
rz(-0.80414786) q[3];
sx q[3];
rz(2.647561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5321396) q[2];
sx q[2];
rz(-2.1707462) q[2];
sx q[2];
rz(-0.54507059) q[2];
rz(2.3405781) q[3];
sx q[3];
rz(-0.89371926) q[3];
sx q[3];
rz(-0.092122294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545749) q[0];
sx q[0];
rz(-1.6681404) q[0];
sx q[0];
rz(0.45904485) q[0];
rz(1.1162988) q[1];
sx q[1];
rz(-0.79367677) q[1];
sx q[1];
rz(2.3150516) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3497884) q[0];
sx q[0];
rz(-1.1513745) q[0];
sx q[0];
rz(-0.38695199) q[0];
x q[1];
rz(0.097657875) q[2];
sx q[2];
rz(-1.6689166) q[2];
sx q[2];
rz(-0.94933921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7956808) q[1];
sx q[1];
rz(-1.7670146) q[1];
sx q[1];
rz(1.0190585) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8070368) q[3];
sx q[3];
rz(-1.6826731) q[3];
sx q[3];
rz(2.3528683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7894342) q[2];
sx q[2];
rz(-2.9106079) q[2];
sx q[2];
rz(-2.3918772) q[2];
rz(-2.0189144) q[3];
sx q[3];
rz(-1.9176982) q[3];
sx q[3];
rz(-2.1813755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4753251) q[0];
sx q[0];
rz(-0.98639494) q[0];
sx q[0];
rz(3.0295897) q[0];
rz(0.22459596) q[1];
sx q[1];
rz(-1.1677531) q[1];
sx q[1];
rz(-0.78132838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3488779) q[0];
sx q[0];
rz(-1.166806) q[0];
sx q[0];
rz(-3.1125665) q[0];
x q[1];
rz(1.4766598) q[2];
sx q[2];
rz(-1.899666) q[2];
sx q[2];
rz(0.48966792) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8857552) q[1];
sx q[1];
rz(-1.9966049) q[1];
sx q[1];
rz(-2.366723) q[1];
rz(-0.16181337) q[3];
sx q[3];
rz(-2.1573632) q[3];
sx q[3];
rz(-0.88255054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3096699) q[2];
sx q[2];
rz(-2.2396542) q[2];
sx q[2];
rz(-0.93264467) q[2];
rz(1.0798838) q[3];
sx q[3];
rz(-2.6974758) q[3];
sx q[3];
rz(3.1058969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2340853) q[0];
sx q[0];
rz(-1.1312753) q[0];
sx q[0];
rz(1.3908516) q[0];
rz(2.8098409) q[1];
sx q[1];
rz(-1.876588) q[1];
sx q[1];
rz(-1.5914241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0858147) q[0];
sx q[0];
rz(-1.5356488) q[0];
sx q[0];
rz(-0.9279535) q[0];
rz(-1.6076902) q[2];
sx q[2];
rz(-2.0301182) q[2];
sx q[2];
rz(1.5806035) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4188288) q[1];
sx q[1];
rz(-2.9420442) q[1];
sx q[1];
rz(-2.2390689) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.329439) q[3];
sx q[3];
rz(-0.5049754) q[3];
sx q[3];
rz(-0.92588378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8288237) q[2];
sx q[2];
rz(-0.88461107) q[2];
sx q[2];
rz(-1.5213607) q[2];
rz(2.17365) q[3];
sx q[3];
rz(-1.5827551) q[3];
sx q[3];
rz(0.91606417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5172326) q[0];
sx q[0];
rz(-2.7150798) q[0];
sx q[0];
rz(2.970001) q[0];
rz(-1.0393556) q[1];
sx q[1];
rz(-1.2587222) q[1];
sx q[1];
rz(-1.4412057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2644669) q[0];
sx q[0];
rz(-1.4466337) q[0];
sx q[0];
rz(-0.32372253) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4588548) q[2];
sx q[2];
rz(-2.0840692) q[2];
sx q[2];
rz(1.7626732) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4413195) q[1];
sx q[1];
rz(-1.1179578) q[1];
sx q[1];
rz(-2.7445566) q[1];
rz(-pi) q[2];
rz(3.0002241) q[3];
sx q[3];
rz(-2.1297567) q[3];
sx q[3];
rz(-1.5822922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8695996) q[2];
sx q[2];
rz(-1.4078434) q[2];
sx q[2];
rz(2.5970411) q[2];
rz(-2.6416687) q[3];
sx q[3];
rz(-2.1130424) q[3];
sx q[3];
rz(-2.669529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8498103) q[0];
sx q[0];
rz(-2.6047459) q[0];
sx q[0];
rz(-2.612402) q[0];
rz(0.57580194) q[1];
sx q[1];
rz(-1.878783) q[1];
sx q[1];
rz(-2.0226488) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7125583) q[0];
sx q[0];
rz(-0.055744113) q[0];
sx q[0];
rz(-0.38614614) q[0];
rz(2.7706657) q[2];
sx q[2];
rz(-1.7417522) q[2];
sx q[2];
rz(1.2990739) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3333867) q[1];
sx q[1];
rz(-1.1115326) q[1];
sx q[1];
rz(0.65826584) q[1];
x q[2];
rz(-1.0435853) q[3];
sx q[3];
rz(-1.7526748) q[3];
sx q[3];
rz(-1.393569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90073663) q[2];
sx q[2];
rz(-2.6764328) q[2];
sx q[2];
rz(-2.4294803) q[2];
rz(1.6093048) q[3];
sx q[3];
rz(-2.1963547) q[3];
sx q[3];
rz(-1.7942662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3466472) q[0];
sx q[0];
rz(-2.0136588) q[0];
sx q[0];
rz(-2.0704863) q[0];
rz(1.2062997) q[1];
sx q[1];
rz(-1.0164398) q[1];
sx q[1];
rz(-1.4869022) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39098719) q[0];
sx q[0];
rz(-1.3840527) q[0];
sx q[0];
rz(-0.67573993) q[0];
rz(-pi) q[1];
rz(-2.5096312) q[2];
sx q[2];
rz(-2.2990531) q[2];
sx q[2];
rz(1.2974844) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8739972) q[1];
sx q[1];
rz(-1.2800299) q[1];
sx q[1];
rz(1.7520755) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28122357) q[3];
sx q[3];
rz(-2.0042002) q[3];
sx q[3];
rz(-3.0126257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6335166) q[2];
sx q[2];
rz(-2.241892) q[2];
sx q[2];
rz(3.0637975) q[2];
rz(0.22732321) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(1.0556861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28854293) q[0];
sx q[0];
rz(-0.74497861) q[0];
sx q[0];
rz(0.37260923) q[0];
rz(-2.695072) q[1];
sx q[1];
rz(-1.6852854) q[1];
sx q[1];
rz(-2.2850697) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73222173) q[0];
sx q[0];
rz(-1.6147095) q[0];
sx q[0];
rz(3.0847021) q[0];
x q[1];
rz(-0.90214731) q[2];
sx q[2];
rz(-2.9716316) q[2];
sx q[2];
rz(2.4580177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9042272) q[1];
sx q[1];
rz(-1.6886535) q[1];
sx q[1];
rz(1.1706283) q[1];
rz(-pi) q[2];
rz(2.4803745) q[3];
sx q[3];
rz(-2.2228129) q[3];
sx q[3];
rz(0.87424247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38241688) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(0.55553931) q[2];
rz(-2.7555452) q[3];
sx q[3];
rz(-1.6470393) q[3];
sx q[3];
rz(-2.4773795) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1356708) q[0];
sx q[0];
rz(-1.6914524) q[0];
sx q[0];
rz(1.2941262) q[0];
rz(2.338943) q[1];
sx q[1];
rz(-1.6812656) q[1];
sx q[1];
rz(0.38801286) q[1];
rz(2.8254208) q[2];
sx q[2];
rz(-2.7718622) q[2];
sx q[2];
rz(-2.3864701) q[2];
rz(-0.070128154) q[3];
sx q[3];
rz(-0.84134103) q[3];
sx q[3];
rz(0.76920912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
