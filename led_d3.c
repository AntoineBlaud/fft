// source /opt/poky-atmel/1.7.1/environment-setup-cortexa5t2hf-vfp-neon-poky-linux-gnueabi

// $CC led_blink.c -o led_blink

// sudo minicom -D /dev/ttyACM0 -b 115200

// mount /dev/sda1 /media

#include <stdio.h> // pour printf

#include <stdlib.h> // pour exit

#include <fcntl.h> // pour open, read, write et close

#include <unistd.h> // pour sleep

#include <errno.h> // pour perror

#define LED_PATH "/sys/class/leds/d3/brightness"

int main(void)
{

    int led_file, i;

    if ((led_file = open(LED_PATH, O_RDWR)) == -1)
    {

        perror(LED_PATH);

        exit(EXIT_FAILURE);
    }
    int mod, led = 0;
    // valeurs possible de led :   '0' -> éteint
    //                             '1' -> allumé
    while (1)
    {
        mod++;
        led = mod % 2 + '0';
        write(led_file, &led, sizeof(char));
        sleep(1);
    }
    close(led_file);

    return EXIT_SUCCESS;
}