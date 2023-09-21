OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(0.8843511) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(1.6391099) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76915574) q[0];
sx q[0];
rz(-0.40524846) q[0];
sx q[0];
rz(-2.0408819) q[0];
rz(-pi) q[1];
rz(-2.5764478) q[2];
sx q[2];
rz(-1.5663212) q[2];
sx q[2];
rz(-0.43484136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.68617649) q[1];
sx q[1];
rz(-1.505291) q[1];
sx q[1];
rz(2.8951365) q[1];
x q[2];
rz(-2.2803335) q[3];
sx q[3];
rz(-1.2951295) q[3];
sx q[3];
rz(2.0303126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.78645906) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(-2.4856429) q[2];
rz(1.9338699) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(-0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(2.7080652) q[0];
rz(0.22878376) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(-3.1343592) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9911306) q[0];
sx q[0];
rz(-0.40369243) q[0];
sx q[0];
rz(1.3315721) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8859293) q[2];
sx q[2];
rz(-1.8506146) q[2];
sx q[2];
rz(-1.8880106) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0011598) q[1];
sx q[1];
rz(-1.3507441) q[1];
sx q[1];
rz(-0.37186719) q[1];
rz(2.1528483) q[3];
sx q[3];
rz(-2.3855004) q[3];
sx q[3];
rz(-2.3088629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6146415) q[2];
sx q[2];
rz(-2.3336637) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(-2.8377623) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4678629) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(1.7720222) q[0];
rz(1.2415775) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(0.26161584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3900671) q[0];
sx q[0];
rz(-1.5528423) q[0];
sx q[0];
rz(-3.1319588) q[0];
rz(-pi) q[1];
rz(-1.4357655) q[2];
sx q[2];
rz(-2.3927852) q[2];
sx q[2];
rz(2.7111862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1640472) q[1];
sx q[1];
rz(-0.79037145) q[1];
sx q[1];
rz(-0.17069823) q[1];
rz(-pi) q[2];
rz(1.5393799) q[3];
sx q[3];
rz(-1.0580225) q[3];
sx q[3];
rz(0.98785066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6802784) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(-0.20351163) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(0.27954277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.571636) q[0];
rz(1.0034026) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(1.8932231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6471449) q[0];
sx q[0];
rz(-1.5054504) q[0];
sx q[0];
rz(1.6676184) q[0];
x q[1];
rz(-2.4291971) q[2];
sx q[2];
rz(-0.27370307) q[2];
sx q[2];
rz(-1.7143539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6971671) q[1];
sx q[1];
rz(-2.4788692) q[1];
sx q[1];
rz(-0.18130937) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8321886) q[3];
sx q[3];
rz(-0.86463118) q[3];
sx q[3];
rz(0.9725001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0288329) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(-1.933243) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(-0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355302) q[0];
sx q[0];
rz(-1.0536138) q[0];
sx q[0];
rz(0.77520448) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(-0.90243375) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31061253) q[0];
sx q[0];
rz(-1.5251625) q[0];
sx q[0];
rz(2.5020585) q[0];
x q[1];
rz(-2.6871704) q[2];
sx q[2];
rz(-1.3497958) q[2];
sx q[2];
rz(-2.1511252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.202995) q[1];
sx q[1];
rz(-2.1122167) q[1];
sx q[1];
rz(0.90403892) q[1];
rz(-1.8875214) q[3];
sx q[3];
rz(-1.5734908) q[3];
sx q[3];
rz(0.56841422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9465785) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(-1.9449332) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(-0.50317558) q[0];
rz(-1.4563837) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(-2.9398289) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87529463) q[0];
sx q[0];
rz(-1.504997) q[0];
sx q[0];
rz(-1.6732897) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1226419) q[2];
sx q[2];
rz(-2.0070224) q[2];
sx q[2];
rz(0.052554616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4469874) q[1];
sx q[1];
rz(-1.1519377) q[1];
sx q[1];
rz(2.6161731) q[1];
rz(0.50815177) q[3];
sx q[3];
rz(-2.3414632) q[3];
sx q[3];
rz(-1.2451764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(0.8231419) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(0.87583035) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(3.0840432) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(-0.94271359) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7563815) q[0];
sx q[0];
rz(-0.74059534) q[0];
sx q[0];
rz(-2.5368607) q[0];
rz(0.96037453) q[2];
sx q[2];
rz(-2.865961) q[2];
sx q[2];
rz(-2.9396217) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49589866) q[1];
sx q[1];
rz(-1.8897448) q[1];
sx q[1];
rz(-1.1369399) q[1];
rz(-pi) q[2];
rz(-1.058217) q[3];
sx q[3];
rz(-1.7985385) q[3];
sx q[3];
rz(2.3678569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1192347) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(2.7589202) q[2];
rz(1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(-2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4246178) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(-0.27012816) q[0];
rz(-0.62942901) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2487508) q[0];
sx q[0];
rz(-2.1584956) q[0];
sx q[0];
rz(2.6177004) q[0];
rz(-1.7173041) q[2];
sx q[2];
rz(-0.97587817) q[2];
sx q[2];
rz(2.811424) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.24208454) q[1];
sx q[1];
rz(-0.71166066) q[1];
sx q[1];
rz(2.3942024) q[1];
rz(0.65981221) q[3];
sx q[3];
rz(-1.5471349) q[3];
sx q[3];
rz(0.82994474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55390629) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(-1.2109057) q[2];
rz(2.8816913) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(-0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0652086) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(-0.22924766) q[0];
rz(-2.8385838) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(-1.680826) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0871353) q[0];
sx q[0];
rz(-1.5883049) q[0];
sx q[0];
rz(2.8580335) q[0];
x q[1];
rz(1.3912348) q[2];
sx q[2];
rz(-1.944724) q[2];
sx q[2];
rz(-2.0342846) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8252392) q[1];
sx q[1];
rz(-1.5484527) q[1];
sx q[1];
rz(-0.52492001) q[1];
rz(-0.78482307) q[3];
sx q[3];
rz(-1.2900969) q[3];
sx q[3];
rz(2.7626038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4298657) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(-0.6955859) q[2];
rz(-0.43186489) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(2.4263884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(0.33690548) q[0];
rz(-0.20740549) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(-0.38063231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779553) q[0];
sx q[0];
rz(-1.5008238) q[0];
sx q[0];
rz(-1.3689343) q[0];
rz(-pi) q[1];
rz(0.74110465) q[2];
sx q[2];
rz(-2.0217102) q[2];
sx q[2];
rz(1.8795183) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14643529) q[1];
sx q[1];
rz(-2.6019115) q[1];
sx q[1];
rz(1.9567009) q[1];
rz(-pi) q[2];
rz(-0.59488876) q[3];
sx q[3];
rz(-0.66685646) q[3];
sx q[3];
rz(-2.8951333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0011065817) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(-2.005119) q[2];
rz(3.100637) q[3];
sx q[3];
rz(-0.80871964) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.3363591) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(1.0271172) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(2.3958191) q[2];
sx q[2];
rz(-2.8095828) q[2];
sx q[2];
rz(-2.7381894) q[2];
rz(0.20053486) q[3];
sx q[3];
rz(-2.0024828) q[3];
sx q[3];
rz(1.9867292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
