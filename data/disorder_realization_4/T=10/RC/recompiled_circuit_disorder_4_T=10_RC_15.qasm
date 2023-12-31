OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.714158) q[0];
sx q[0];
rz(-2.7058869) q[0];
sx q[0];
rz(-0.92619196) q[0];
rz(-1.1821094) q[1];
sx q[1];
rz(3.8745772) q[1];
sx q[1];
rz(12.193845) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6719138) q[0];
sx q[0];
rz(-0.9979453) q[0];
sx q[0];
rz(1.3290622) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9973642) q[2];
sx q[2];
rz(-0.63604522) q[2];
sx q[2];
rz(-2.6682105) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2791427) q[1];
sx q[1];
rz(-1.3819873) q[1];
sx q[1];
rz(-0.11629176) q[1];
rz(1.7872693) q[3];
sx q[3];
rz(-1.4604005) q[3];
sx q[3];
rz(2.20761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42840502) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(0.92450809) q[2];
rz(1.472578) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(-1.6424461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9532042) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(1.6054608) q[0];
rz(0.19451441) q[1];
sx q[1];
rz(-1.8201927) q[1];
sx q[1];
rz(-0.054873437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7354436) q[0];
sx q[0];
rz(-2.9950954) q[0];
sx q[0];
rz(-0.88645257) q[0];
rz(-pi) q[1];
rz(0.61972159) q[2];
sx q[2];
rz(-2.0267068) q[2];
sx q[2];
rz(2.1330619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3000852) q[1];
sx q[1];
rz(-1.9133139) q[1];
sx q[1];
rz(-0.0042580982) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67316405) q[3];
sx q[3];
rz(-1.4041956) q[3];
sx q[3];
rz(1.8002321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43869552) q[2];
sx q[2];
rz(-0.56151152) q[2];
sx q[2];
rz(-1.4820209) q[2];
rz(2.1510018) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.7747442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7822587) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(-2.7666132) q[0];
rz(-0.24770501) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.8050271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35555392) q[0];
sx q[0];
rz(-0.22909129) q[0];
sx q[0];
rz(1.8285719) q[0];
rz(-pi) q[1];
rz(-1.1281563) q[2];
sx q[2];
rz(-1.0269564) q[2];
sx q[2];
rz(1.3904872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3805483) q[1];
sx q[1];
rz(-2.3420463) q[1];
sx q[1];
rz(1.9995081) q[1];
rz(-pi) q[2];
rz(-2.7258354) q[3];
sx q[3];
rz(-2.3299179) q[3];
sx q[3];
rz(-0.31879253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(-1.0245727) q[2];
rz(1.2767977) q[3];
sx q[3];
rz(-1.5494616) q[3];
sx q[3];
rz(1.6259441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.9643726) q[0];
sx q[0];
rz(-1.465088) q[0];
sx q[0];
rz(-0.033551034) q[0];
rz(-1.230348) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(-1.7376602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0516292) q[0];
sx q[0];
rz(-1.7647867) q[0];
sx q[0];
rz(-0.35332638) q[0];
rz(-2.3061182) q[2];
sx q[2];
rz(-0.88047853) q[2];
sx q[2];
rz(0.26963216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0714604) q[1];
sx q[1];
rz(-1.6825819) q[1];
sx q[1];
rz(1.2791355) q[1];
x q[2];
rz(-2.2622044) q[3];
sx q[3];
rz(-1.6367568) q[3];
sx q[3];
rz(1.1018167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6491062) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(-0.21326324) q[2];
rz(0.02031859) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7810818) q[0];
sx q[0];
rz(-2.0251944) q[0];
sx q[0];
rz(-1.0158585) q[0];
rz(1.6332743) q[1];
sx q[1];
rz(-2.1834178) q[1];
sx q[1];
rz(0.034084592) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3344454) q[0];
sx q[0];
rz(-2.7140518) q[0];
sx q[0];
rz(-0.39512511) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5675315) q[2];
sx q[2];
rz(-1.3458369) q[2];
sx q[2];
rz(-3.0846734) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2451671) q[1];
sx q[1];
rz(-2.8571269) q[1];
sx q[1];
rz(-1.5462272) q[1];
rz(-0.91655101) q[3];
sx q[3];
rz(-0.74186462) q[3];
sx q[3];
rz(1.786701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4013227) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(1.9704698) q[2];
rz(1.8088388) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(-0.12935054) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4532582) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(-2.7057498) q[0];
rz(2.0360937) q[1];
sx q[1];
rz(-1.8445797) q[1];
sx q[1];
rz(-1.2058535) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57731956) q[0];
sx q[0];
rz(-0.84007971) q[0];
sx q[0];
rz(-2.8217836) q[0];
rz(-1.7370991) q[2];
sx q[2];
rz(-0.096509343) q[2];
sx q[2];
rz(-1.4788747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.3881512) q[1];
sx q[1];
rz(-2.6302164) q[1];
sx q[1];
rz(-3.0877153) q[1];
rz(-pi) q[2];
rz(1.4924963) q[3];
sx q[3];
rz(-2.7408528) q[3];
sx q[3];
rz(2.0164255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2445406) q[2];
sx q[2];
rz(-2.5138833) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(2.7339325) q[3];
sx q[3];
rz(-1.9577273) q[3];
sx q[3];
rz(-0.44529644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62149858) q[0];
sx q[0];
rz(-0.61722732) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(-0.78397059) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(2.6046682) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85978956) q[0];
sx q[0];
rz(-1.4872695) q[0];
sx q[0];
rz(-0.066017166) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0579254) q[2];
sx q[2];
rz(-0.29209902) q[2];
sx q[2];
rz(-2.3592754) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9432224) q[1];
sx q[1];
rz(-2.3343625) q[1];
sx q[1];
rz(-0.13065773) q[1];
rz(-1.4173996) q[3];
sx q[3];
rz(-0.985257) q[3];
sx q[3];
rz(0.034754001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9817104) q[2];
sx q[2];
rz(-1.2572224) q[2];
sx q[2];
rz(-2.9505777) q[2];
rz(-0.28132176) q[3];
sx q[3];
rz(-1.1912991) q[3];
sx q[3];
rz(-2.7991926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1087588) q[0];
sx q[0];
rz(-0.37029752) q[0];
sx q[0];
rz(-0.38761815) q[0];
rz(3.0265813) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(0.33755916) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96466366) q[0];
sx q[0];
rz(-1.3725855) q[0];
sx q[0];
rz(-0.49074178) q[0];
rz(2.2099724) q[2];
sx q[2];
rz(-1.6281623) q[2];
sx q[2];
rz(0.66040874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2355289) q[1];
sx q[1];
rz(-1.7019094) q[1];
sx q[1];
rz(0.99793418) q[1];
x q[2];
rz(-2.3861305) q[3];
sx q[3];
rz(-1.6172234) q[3];
sx q[3];
rz(1.4481627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9541624) q[2];
sx q[2];
rz(-2.607441) q[2];
sx q[2];
rz(1.2672651) q[2];
rz(-2.3665442) q[3];
sx q[3];
rz(-1.6036443) q[3];
sx q[3];
rz(-2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545814) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(2.9206081) q[0];
rz(0.9221319) q[1];
sx q[1];
rz(-1.2698959) q[1];
sx q[1];
rz(2.1386713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5054277) q[0];
sx q[0];
rz(-1.3190077) q[0];
sx q[0];
rz(1.7136784) q[0];
rz(-2.1557751) q[2];
sx q[2];
rz(-2.0009811) q[2];
sx q[2];
rz(0.9790203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50833118) q[1];
sx q[1];
rz(-2.0610626) q[1];
sx q[1];
rz(-1.0432748) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8691611) q[3];
sx q[3];
rz(-2.2386132) q[3];
sx q[3];
rz(-1.2215134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6691436) q[2];
sx q[2];
rz(-0.74206918) q[2];
sx q[2];
rz(-2.9830902) q[2];
rz(-2.6878099) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(0.84428549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52699387) q[0];
sx q[0];
rz(-2.232593) q[0];
sx q[0];
rz(2.9504543) q[0];
rz(0.29516164) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(-0.89231649) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44360456) q[0];
sx q[0];
rz(-2.3058878) q[0];
sx q[0];
rz(-2.9025335) q[0];
rz(-pi) q[1];
rz(-3.0478165) q[2];
sx q[2];
rz(-2.2095592) q[2];
sx q[2];
rz(0.057657777) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3117864) q[1];
sx q[1];
rz(-1.808842) q[1];
sx q[1];
rz(1.6386375) q[1];
x q[2];
rz(1.8252556) q[3];
sx q[3];
rz(-1.99031) q[3];
sx q[3];
rz(-2.6655243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7717379) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(-0.85956335) q[2];
rz(1.2236979) q[3];
sx q[3];
rz(-2.2151291) q[3];
sx q[3];
rz(-2.3971476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201465) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(-1.4795115) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(1.2791469) q[2];
sx q[2];
rz(-1.5970061) q[2];
sx q[2];
rz(2.9562052) q[2];
rz(-1.6395232) q[3];
sx q[3];
rz(-1.5309661) q[3];
sx q[3];
rz(-1.0504709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
